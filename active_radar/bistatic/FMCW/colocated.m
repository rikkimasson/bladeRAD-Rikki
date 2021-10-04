clear all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 15;    % Expeiment Name
capture_duration = 30;        % capture duration
Fs = 40e6;          % Sample Rate of SDR per I & Q (in reality Fs is double this)
pulse_duration = 1e-3;   % Desired Pulse Duration 
Bw = 40e6;          % LFM Bandwidth 
save_directory = "/home/piers/Documents/Captures/"; % rach experiment will save as a new folder in this directory

% Radar Parameters 
Fc = 5800e6;   % Central RF 
Tx_gain = 66;  % [-23.75, 66]
Rx1_gain = 36;
Rx2_gain = 0;
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
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   RF_freq,...
                                   Bw_M,...
                                   Tx_SDR,...
                                   'slave',...
                                   3,...
                                   'tx');
%     tx_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(tx_command);
    pause(5);

    % Setup mono Rx SDR 
    [trig_flag_2,rx_command] = create_shell_command(Experiment_ID,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   RF_freq,...
                                   Bw_M,...
                                   Bi_Rx_SDR,...
                                   'slave',...
                                   1,...
                                   'rx'); 
    if trig_flag_1 && trig_flag_2
        "Trigger Conflict - FMCW Radar"
        return
    end                                                          
    system(rx_command); % Blocking system command execution 

    
    % Setup Pseudo Bistatic Rx SDR 
    [trig_flag_3,rx_command] = create_shell_command(Experiment_ID,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   RF_freq,...
                                   Bw_M,...
                                   Rx_SDR,...
                                   'master',...
                                   1,...
                                   'rx'); 
    if trig_flag_1 && trig_flag_2 && trig_flag_3
        "Trigger Conflict - FMCW Radar"
        return
    end                                                          
    system(rx_command); % Blocking system command execution     

%% Save Raw Data and create  header to directory 
    exp_dir = save_directory + Experiment_ID + '/';
    make_dir = 'mkdir ' + exp_dir;
    system(make_dir); % Blocking system command execution
    move_file = 'mv /tmp/fmcw_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Rx Data Copyied to Save directory"
    else 
        "Rx Copy Failed"
        return
        
    end
    save(exp_dir + 'FMCW Experimental Configuration') 

    
    
%% Load Reference Deramp Signal
    refsig = load_refsig(Bw_M,Fc,pulse_duration);
    %figure 
    %spectrogram(refsig,128,100,100,Fs,'centered','yaxis')

    
    
%% Load Signal, Mix and Dermap Signal  
file_location = exp_dir + 'fmcw_' + Experiment_ID;
[max_range_actual,processed_signal] = deramp_and_decimate(file_location,max_range,refsig,capture_duration,number_pulses,Fs,slope);
save(exp_dir + 'deramped_signal','processed_signal')


%% Plot RTI

    Range_axis = linspace(0,max_range_actual,size(processed_signal,1));
    Range_bin = 1:size(processed_signal,1);
    time_axis = linspace(0,size(processed_signal,2)*pulse_duration,size(processed_signal,2));
    RTI_plot= transpose(10*log10(abs(processed_signal./max(processed_signal(:)))));
    figure
    imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
    % % % xlim([0 25])
    % % %ylim([0 0.0005])
    % % grid on            
    % % colorbar
    % % ylabel('Time (Sec)')
    % % xlabel('Range (m)')   
    % % fig_title = "Monostatic RTI - Test " + Test_id;
    % % title(fig_title);fig_name = save_directory + "/RTI_" + Test_id + ".jpg";
    % % saveas(fig,fig_name,'jpeg') 
    % % % fig_name = "/home/piers/Desktop/FMCW/Experiments/Test_"+ Test_id + ".fig";
    % % % savefig(fig_name)

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