clear all
close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 1;       % Expeiment Name
capture_duration = 5;    % capture duration
Fs = 30e6;               % Sample Rate of SDR per I & Q (in reality Fs is double this)
pulse_duration = 0.5e-3;   % Desired Pulse Duration 
Bw = 30e6;               % LFM Bandwidth 
%save_directory = "/media/sdrlaptop1/T7/22_06_21_N0/"; % each experiment will save as a new folder in this directory
save_directory = "~/Documents/bladeRAD_Captures/06_07_2022_farm/multistatic/s_band/"; % each experiment will save as a new folder in this directory
exp_dir = save_directory + Experiment_ID + '/';

% Radar Parameters 
Fc = 2.44e9;   % Central RF 
Tx_gain = 60;  % [-23.75, 66] (S-Band = 23.5 dBm) (C-Band = 15.8 dBm)
Rx1_gain = 0;  % [-16, 60]
Rx2_gain = 0;  % [-16, 60]
Tx_SDR = 1;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
Rx_SDR = 2;   % SDR to use for RX

% Procesing Parameters
 max_range = 1000; %max range to LPF filter data to
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
chirp = saw_LFM_chirp(Bw,pulse_duration,Fs);
% chirp = tone_generator(-10e6,pulse_duration,Fs);
% 
% figure
% plot(real(chirp));
% spectrogram(chirp)
% spec= (fft(chirp));
% figure
% plot(linspace(-15e6,15e6,length(spec)),20*log10(abs(spec)))
save_sc16q11('/tmp/chirp.sc16q11', chirp); %save chirp to binary file
% clear chirp
    % spectrogram(chirp,128,100,128,Fs,'centered','yaxis') %plot spectrogram of chirp
    % f = linspace(-0.5 * Fs, 0.5 * Fs, length(chirp));
    % figure
    % plot(f,20*log10(abs(fftshift(fft(chirp)))/Samps_per_pulse));%plot FFT of chirp to show entire spectral content

%% Setup Radar
    % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''

    % Setup Tx SDR 
    [trig_flag_1,tx_command] = create_shell_command( ...
                                   true,...
                                   Experiment_ID,...
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
    [trig_flag_2,rx_command] = create_shell_command( ...
                                   false,...
                                   Experiment_ID,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   0,...                                 
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   Fc_M,...
                                   Bw_M,...
                                   Rx_SDR,...
                                   'master',...
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
    refsig = load_refsig(Bw_M,Fs,Fc,pulse_duration);
%     figure 
%     spectrogram(refsig,128,100,100,Fs,'centered','yaxis')
    
%% Load Signal, Mix and Dermap Signal  
zero_padding = 1;
lp_filter = false;
file_location = exp_dir + 'active_' + Experiment_ID;
[max_range_actual,deramped_signal] = deramp_and_decimate(file_location,max_range,refsig,capture_duration,number_pulses,Fs,slope,lp_filter);
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
    RTI_plot= transpose(20*log10(abs(processed_signal./max(processed_signal(:)))));
    figure
    fig = imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
%         xlim([0 100])
        grid on            
        colorbar
        ylabel('Time (Sec)')
        xlabel('Range (m)')   
        title("FMCW RTI - " + Experiment_ID);
        fig_name = exp_dir + "FMCW_RTI_" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg') 
% % MTI RTI
%     Range_axis = linspace(0,max_range_actual,size(processed_signal,1));
%     Range_bin = 1:size(processed_signal,1);
%     time_axis = linspace(0,size(processed_signal,2)*pulse_duration,size(processed_signal,2));
%     RTI_plot= transpose(20*log10(abs(MTI_Data./max(MTI_Data(:)))));
%     figure
%     fig = imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
%         xlim([0 100])
%         grid on            
%         colorbar
%         ylabel('Time (Sec)')
%         xlabel('Range (m)')   
%         title("MTI FMCW RTI - " + Experiment_ID);
%         fig_name = exp_dir + "MTI_FMCW_RTI_" + Experiment_ID + ".jpg";
%         saveas(fig,fig_name,'jpeg') 

% plot a single pulse
    figure
    fig = plot(Range_axis,RTI_plot(10,:));
        title("Single Pulse - " + Experiment_ID);
%         xlim([0 100])
        grid on
        ylabel('Relative Power (dB)')
        xlabel('Range (m)')  
        fig_name = exp_dir + "Single_Pulse" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg') 

 % plot phase of maximum return
    [value, pulse_rbin] = max(processed_signal(:,100)); 
    phase_8 = angle(processed_signal(pulse_rbin,:));
    phase_norm_8 = phase_8 - mean(phase_8(1:100));
    phase_time = phase_8/(2*pi*Fc);
    figure
    fig = plot(time_axis,phase_8)
    xlabel('Time (s)')
    ylabel('phase (radians)')
    fig_name = exp_dir + "Phase Series" + Experiment_ID + ".jpg";
    saveas(fig,fig_name,'jpeg') 
    
    nfft = 1024;
    figure
    f = linspace(-0.5 * 1000, 0.5 * 1000, nfft);
    freq_dom = fftshift(abs(fft(processed_signal(pulse_rbin,:),nfft)));
    plot(f,20*log10(freq_dom./max(freq_dom(:))))
    grid on; grid minor;
    xlabel('Frequency');
    ylabel('Relative Power (dB)')  
    
    hold on;
    [pow] = pwelch(processed_signal(pulse_rbin,:),hamming(nfft),0.9,'centered','power');
%     figure    
    f = linspace(-0.5 * 2000, 0.5 * 2000, nfft);
    plot(f,10*log10(pow./max(pow(:))))
    grid on; grid minor;
    xlabel('Frequency (Hz)');
    ylabel('Relative Power (dB)')      
             
%% Plot Spectrogram 
    r_start = 1;
    r_stop = 100;
    l_fft = 1024;
    pad_factor = 1;
    overlap_factor = 0.99;
    integrated_data = sum(processed_signal(r_start:r_stop,:));

    [spect,f] = spectrogram(integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-50 0]);
        ylim([-100 100])
        colorbar
        xlabel('Time (Sec)')
        % ylabel('Radial Velocity (mph)')   
        ylabel('Doppler Frequency (Hz)')  
        fig_title = "Spectrogram - Test " + Experiment_ID;
        title(fig_title);
        fig_name = exp_dir + "/Spectrogram_" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')        

  %Plot MTI Spectrogram 
    mti_integrated_data = sum(MTI_Data(r_start:r_stop,:));
    [spect,f] = spectrogram(mti_integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-50 0]);
        ylim([-100 100])
        colorbar
        xlabel('Time (Sec)')
        % ylabel('Radial Velocity (mph)')   
        ylabel('Doppler Frequency (Hz)')  
        fig_title = "Monostatic Single Delay Line MTI Spectrogram - Test " + Experiment_ID;
        title(fig_title);
        fig_name = exp_dir + "/MTI_Spectrogram_" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')        