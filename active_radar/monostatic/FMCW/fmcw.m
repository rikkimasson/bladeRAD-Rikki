clear all
close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 9;       % Expeiment Name
capture_duration = 1;    % capture duration
Fs = 20e6;               % Sample Rate of SDR per I & Q (in reality Fs is double this)
pulse_duration = 1e-3;   % Desired Pulse Duration 
Bw = 20e6;               % LFM Bandwidth 
save_directory = "/media/sdrlaptop1/T7/22_06_21_N0/"; % each experiment will save as a new folder in this directory
    exp_dir = save_directory + Experiment_ID + '/';

% Radar Parameters 
Fc = 2.4e9;   % Central RF 
Tx_gain = 30;  % [-23.75, 66] (S-Band = 23.5 dBm) (C-Band = 15.8 dBm)
Rx1_gain = 26;  % [-16, 60]
Rx2_gain = 0;  % [-16, 60]
Tx_SDR = 1;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
Rx_SDR = 2;   % SDR to use for RX

% Procesing Parameters
 max_range = 100; %max range to LPF filter data to
 processing_flag = true;

% Parameters not configurable by user 
    C = physconst('LightSpeed');
    PRF = 1/pulse_duration;
    slope = Bw/pulse_duration;
    F_Max = Bw/2;
    R_Max = beat2range(F_Max,slope);
    sample_duration = 1/Fs;
    samples_per_pulse = pulse_duration/sample_duration
    number_pulses = capture_duration / pulse_duration
    number_cap_samps = capture_duration/sample_duration
    Fc_M = Fc/1e6;   % RF in MHz 
    Bw_M = Bw/1e6;      % BW in MHz


%% Create Sawtooth Chirp for bladeRF
refsig = saw_LFM_chirp(Bw,pulse_duration,Fs);
figure
plot(real(chirp));
save_sc16q11('/tmp/chirp.sc16q11', chirp); %save chirp to binary file
clear chirp
    spectrogram(chirp,128,100,128,Fs,'centered','yaxis') %plot spectrogram of chirp
    f = linspace(-0.5 * Fs, 0.5 * Fs, length(chirp));
    figure
    plot(f,20*log10(abs(fftshift(fft(chirp)))/size(chirp,1)));%plot FFT of chirp to show entire spectral content

%% Setup Radar
    % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''

    % Setup Tx SDR 
    [trig_flag_1,tx_command] = create_shell_command(Experiment_ID,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   0,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   Fc_M,...
                                   Bw_M,...
                                   Tx_SDR,...
                                   'slave',...
                                   3,...
                                   'tx');
    tx_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(tx_command);
    pause(5);

    % Setup Rx SDR 
    [trig_flag_2,rx_command] = create_shell_command(Experiment_ID,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   0,...                                 
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   Fc_M,...
                                   Bw_M,...
                                   Rx_SDR,...
                                   'slave',...
                                   1,...
                                   'rx'); 
    if trig_flag_1 && trig_flag_2
        "Trigger Conflict - FMCW Radar"
        return
    end                                                          
    system(rx_command) % Blocking system command execution 


%% Save Raw Data and create  header to directory 
    make_dir = 'mkdir ' + exp_dir;
    system(make_dir); % Blocking system command execution
    move_file = 'mv /tmp/active_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Rx Data Copyied to Save directory"
    else 
        "Rx Copy Failed"
        return
        
    end
    save(exp_dir + 'Experimental Configuration') 

if processing_flag == false
    return 
end
    
%% Load Reference Deramp Signal
    refsig = load_refsig(Bw_M,Fc,pulse_duration);
%     figure 
%     spectrogram(refsig,128,100,100,Fs,'centered','yaxis')
    
%% Load Signal, Mix and Dermap Signal  
zero_padding = 2;
file_location = exp_dir + 'active_' + Experiment_ID;
[max_range_actual,deramped_signal] = deramp_and_decimate(file_location,max_range,refsig,capture_duration,number_pulses,Fs,slope);
save(exp_dir + 'deramped_signal','deramped_signal')


%% Window and FFT Signal 
    % window signal
        w = window('hann',size(deramped_signal,1));
        windowed_signal = deramped_signal.*w;
    % fft signal
        zero_padding = 1;
        processed_signal = fft(windowed_signal,size(windowed_signal,1)*zero_padding);
    
    
    %% MTI Filtering 
        % Single Delay Line Filter 
        MTI_Data = zeros(size(processed_signal));
              for i=2:number_pulses
                    MTI_Data(:,i) = processed_signal(:,i)-processed_signal(:,i-1);
              end



%% Plot RTI

    Range_axis = linspace(0,max_range_actual,size(processed_signal,1));
    Range_bin = 1:size(processed_signal,1);
    time_axis = linspace(0,size(processed_signal,2)*pulse_duration,size(processed_signal,2));
    RTI_plot= transpose(10*log10(abs(processed_signal./max(processed_signal(:)))));
    figure
    fig = imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
        xlim([0 100])
        grid on            
        colorbar
        ylabel('Time (Sec)')
        xlabel('Range (m)')   
        title("FMCW RTI - " + Experiment_ID);
        fig_name = exp_dir + "FMCW_RTI_" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg') 


    figure
    plot(Range_axis,RTI_plot(10,:));
        title("Single Pulse - " + Experiment_ID);
        xlim([0 100])
        grid on
        ylabel('Relative Power (dB)')
        xlabel('Range (m)')  
        fig_name = exp_dir + "Single_Pulse" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg') 

[value, pulse_rbin] = max(processed_signal(:,100)); 
phase_9 = angle(processed_signal(pulse_rbin,:));
phase_norm_9 = phase_9 - mean(phase_9(1:100));
save(exp_dir + 'phase_norm_9', 'phase_9')
figure
plot(phase_norm_9)
figure
plot(abs(processed_signal(pulse_rbin,:)))
        
        
        %         
% % Plot Spectrogram 
%     r_bin = 3;
%     l_fft = 1024;
%     pad_factor = 1;
%     overlap_factor = 0.99;
%     [spect,f] = spectrogram(processed_signal(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
%         % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
%         v=dop2speed(f,C/Fc)*2.237;
%         spect= 10*log10(abs(spect./max(spect(:))));
%         figure
%         fig = imagesc(time_axis,f,spect,[-50 0]);
%         ylim([-100 100])
%         colorbar
%         xlabel('Time (Sec)')
%         % ylabel('Radial Velocity (mph)')   
%         ylabel('Doppler Frequency (Hz)')  
%         fig_title = "Monostatic Single Delay Line MTI Spectrogram - Test " + Experiment_ID;
%         title(fig_title);
%         fig_name = save_directory + "/MTI_Spectrogram_" + Experiment_ID + ".jpg";
%         saveas(fig,fig_name,'jpeg')        
%         
% %% MTI Filtering 
%     % Single Delay Line Filter 
%     MTI_Data = zeros(size(processed_signal));
%           for i=2:number_pulses
%                 MTI_Data(:,i) = processed_signal(:,i)-processed_signal(:,i-1);
%           end
%       
%     % Plot MTI RTI      
%     MTI_RTI_plot= transpose(10*log10(abs(MTI_Data./max(MTI_Data(:)))));
%     figure
%     fig = imagesc(Range_bin,time_axis,MTI_RTI_plot,[-30,0]);
%         xlim([1 50])
%         %ylim([0 0.0005])
%         grid on            
%         colorbar
%         ylabel('Time (Sec)')
%         xlabel('Range Bin')   
%         fig_title = "Monostatic Single Delay Line MTI  RTI - Test " + Experiment_ID;
%         title(fig_title);
%         fig_name = save_directory + "/MTI_RTI_" + Experiment_ID + ".jpg";
%         saveas(fig,fig_name,'jpeg')
%         plot_signal = toc     
%     
% 
%     %Plot MTI Spectrogram 
%     r_bin = 3;
%     l_fft = 512;
%     pad_factor = 4;
%     overlap_factor = 0.99;
%     [spect,f] = spectrogram(MTI_Data(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
%         % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
%         v=dop2speed(f,C/Fc)*2.237;
%         spect= 10*log10(abs(spect./max(spect(:))));
%         figure
%         fig = imagesc(time_axis,f,spect,[-20 0]);
%         ylim([-100 100])
%         colorbar
%         xlabel('Time (Sec)')
%         % ylabel('Radial Velocity (mph)')   
%         ylabel('Doppler Frequency (Hz)')  
%         fig_title = "Monostatic Single Delay Line MTI Spectrogram - Test " + Experiment_ID;
%         title(fig_title);
%         fig_name = save_directory + "/MTI_Spectrogram_" + Experiment_ID + ".jpg";
%         saveas(fig,fig_name,'jpeg')        

        
%     %% Coherent integration 
%     compressed_data = sum(Dec_Deramped,2);
%     figure
%     plot(Range_axis,abs(compressed_data))
%     xlim([0 1000])

