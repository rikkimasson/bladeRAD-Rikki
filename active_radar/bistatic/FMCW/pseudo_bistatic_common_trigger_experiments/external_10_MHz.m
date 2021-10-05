clear all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 18;    % Expeiment Name
capture_duration = 30;        % capture duration
Fs = 20e6;          % Sample Rate of SDR per I & Q (in reality Fs is double this)
pulse_duration = 1e-3;   % Desired Pulse Duration 
Bw = 20e6;          % LFM Bandwidth 
save_directory = "/home/piers/Documents/Captures/multistatic/"; % rach experiment will save as a new folder in this directory

% Radar Parameters 
Fc = 2400e6;   % Central RF 
Tx_gain = 50  % [-23.75, 66]
Rx1_gain = 36;
Rx2_gain = 36;
Tx_SDR = 1;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
Rx_SDR = 2;   % SDR to use for RX
Bi_Rx_SDR = 3;

% Procesing Parameters
 max_range = 1000; %max range to LPF filter data to

% Parameters not configurable by user 
    C = physconst('LightSpeed');
    PRF = 1/pulse_duration;
    slope = Bw/pulse_duration;
    F_Max = Bw/2;
    R_Max = beat2range(F_Max,slope);
    sample_duration = 1/Fs;
    samples_per_pulse = pulse_duration/sample_duration;
    number_pulses = capture_duration / pulse_duration;
    number_cap_samps = capture_duration/sample_duration;
    RF_freq = Fc/1e6;   % RF in MHz 
    Bw_M = Bw/1e6;      % BW in MHz


%% Create Sawtooth Chirp for bladeRF
chirp = saw_LFM_chirp(Bw,pulse_duration,Fs);
save_sc16q11('/tmp/chirp.sc16q11', chirp); %save chirp to binary file
clear chirp
    % spectrogram(chirp,128,100,128,Fs,'centered','yaxis') %plot spectrogram of chirp
    % f = linspace(-0.5 * Fs, 0.5 * Fs, length(chirp));
    % figure
    % plot(f,20*log10(abs(fftshift(fft(chirp)))/Samps_per_pulse));%plot FFT of chirp to show entire spectral content

%% Setup Radar
    % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''

    % Setup mono Tx SDR 
    [trig_flag_1,tx_command] = create_shell_command(Experiment_ID,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   0,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   RF_freq,...
                                   Bw_M,...
                                   Tx_SDR,...
                                   'slave',...
                                   2,...
                                   'tx');
    tx_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(tx_command);
    pause(5);

    % Setup mono Rx SDR 
    [trig_flag_2,mono_rx_command] = create_shell_command(Experiment_ID,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   0,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   RF_freq,...
                                   Bw_M,...
                                   Rx_SDR,...
                                   'slave',...
                                   1,...
                                   'rx'); 
    if trig_flag_1 && trig_flag_2
        "Trigger Conflict - FMCW Radar"
        return
    end        
    mono_rx_command = mono_rx_command + "&"; % uncomment for non-blocking system command execution
    system(mono_rx_command); % Blocking system command execution 
    pause(5);

    
    % Setup Pseudo Bistatic Rx SDR 
    bi_name = "bistatic_" + Experiment_ID;
    [trig_flag_3,bi_rx_command] = create_shell_command(bi_name,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   0,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   RF_freq,...
                                   Bw_M,...
                                   Bi_Rx_SDR,...
                                   'master',...
                                   2,...
                                   'rx'); 
    if trig_flag_1 && trig_flag_2 && trig_flag_3
        "Trigger Conflict - FMCW Radar"
        return
    end                                                          
    system(bi_rx_command); % Blocking system command execution     

%% Save Raw Data and create  header to directory 
    exp_dir = save_directory + Experiment_ID + '/';
    make_dir = 'mkdir ' + exp_dir;
    system(make_dir); % Blocking system command execution
    move_file = 'mv /tmp/active_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Mono Rx Data Copyied to Save directory"
    else 
        "Mono Rx Copy Failed"
        return  
    end
    % move bi file
    move_file = 'mv /tmp/active_' + string(bi_name) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Bi Rx Data Copyied to Save directory"
    else 
        "Bi Rx Copy Failed"
        return    
    end
    save(exp_dir + 'FMCW Experimental Configuration'); 

    
    
%% Load Reference Deramp Signal
    refsig = load_refsig(Bw_M,Fc,pulse_duration);
    %figure 
    %spectrogram(refsig,128,100,100,Fs,'centered','yaxis')

    
    
%% Load Signal, Mix and Dermap Signal  
zero_padding = 1;
% Monostatic processing 
    file_location = exp_dir + 'active_' + Experiment_ID;
    [max_range_actual,mono_processed_signal] = deramp_and_decimate(file_location,max_range,refsig,capture_duration,number_pulses,Fs,slope,zero_padding);
% Bistatic processing 
    file_location = exp_dir + 'active_' + bi_name;
    [max_range_actual,bi_processed_signal] = deramp_and_decimate(file_location,max_range,refsig,capture_duration,number_pulses,Fs,slope,zero_padding);
    save(exp_dir + 'deramped_signal','mono_processed_signal','bi_processed_signal')
 



%% Plot RTI
    % create axis
        Range_axis = linspace(0,max_range_actual,size(mono_processed_signal,1));
        Range_bin = 1:size(mono_processed_signal,1);
        time_axis = linspace(0,size(mono_processed_signal,2)*pulse_duration,size(mono_processed_signal,2));
    % create mono RTI
        mono_RTI_plot = transpose(10*log10(abs(mono_processed_signal./max(mono_processed_signal(:)))));
    % create bi RTI
        bi_RTI_plot = transpose(10*log10(abs(bi_processed_signal./max(bi_processed_signal(:)))));
    % plot figure
    figure
    fig = subplot(1,2,1);
        imagesc(Range_bin,time_axis,mono_RTI_plot,[-20,0]);   
        xlim([0 100])
        grid on            
        colorbar
        ylabel('Time (Sec)')
        xlabel('Range (m)')   
        fig_title = "Monostatic RTI - Test " + Experiment_ID;
        title(fig_title)
    hold on
    subplot(1,2,2);
        imagesc(Range_bin,time_axis,bi_RTI_plot,[-20,0]);   
        xlim([0 100])
        grid on            
        colorbar
        ylabel('Time (Sec)')
        xlabel('Range (m)')   
        fig_title = "Bistatic RTI - Test " + Experiment_ID;
        title(fig_title)
     fig_name = exp_dir + "RTI_" + Experiment_ID + ".jpg";
     saveas(fig,fig_name,'jpeg')
     fig_name = exp_dir + "RTI_" + Experiment_ID;
     saveas(fig,fig_name)

 % plot phase    
figure
plot(time_axis ,angle(mono_processed_signal(4,:)))
hold on 
fig = plot(time_axis ,angle(bi_processed_signal(4,:)))
        legend("Monostatic Phase","Bistatic Phase")
        xlabel('Time (s)')
        ylabel('Phase (radians)')   
        grid on
        fig_name = exp_dir + "Phase_" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')
        fig_name = exp_dir + "Phase_" + Experiment_ID;
        saveas(fig,fig_name)

%     %% Coherent integration 
%     compressed_data = sum(Dec_Deramped,2);
%     figure
%     plot(Range_axis,abs(compressed_data))
%     xlim([0 1000])
% 
%     %% Spectrogram 
%     spec_bin = 3;
%     l_fft = 512;
%     pad_factor = 4;
%     overlap_factor = 0.99;
%     [spect,f] = spectrogram(Final_Data(spec_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
%     % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
%     v=dop2speed(f,C/Fc)*2.237;
%     spect= 10*log10(abs(spect./max(spect(:))));
%     figure
%     fig = imagesc(time_axis,f,spect,[-50 0]);   
%     ylim([-100 100])
%     colorbar
%     xlabel('Time (Sec)')
%     % ylabel('Radial Velocity (mph)')   
%     ylabel('Doppler Frequency (Hz)')  
%     fig_title = "Monostatic Spectrogram - Test " + Test_id;
%     title(fig_title);
%     fig_name = save_directory + "/Spectrogram_" + Test_id + ".jpg";
%     saveas(fig,fig_name,'jpeg')
