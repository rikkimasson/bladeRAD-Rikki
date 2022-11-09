% clear all
% % close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 1;    % Expeiment Name
capture_duration = 2;        % capture duration
save_directory = "/media/piers/T7/FMCW_Loopback_Range_Calibration"; % each experiment will save as a new folder in this directory
exp_dir = save_directory + Experiment_ID + '/';


% FMCW Parameters 
active.Fs = 30e6;          % Sample Rate of SDR per I & Q (in reality Fs is double this)
active.pulse_duration = 1e-3;   % Desired Pulse Duration 
active.Bw = 30e6;          % LFM Bandwidth 
active.Fc = 2.44e9;   % Central RF 
active.Tx_gain = 60;       % [-23.75, 66] (S-Band = 23.5 dBm) (C-Band = 15.8 dBm)
active.Rx1_gain = 0;      % [-16, 60]
active.Rx2_gain = 0;       % [-16, 60]
Rx_1_lna = true;
Rx_2_lna = true;
active.Tx_SDR = 1;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
active.Rx_SDR = 2;   % SDR to use for RX
% Procesing Parameters
active.max_range = 2000; %max range to LPF filter data to

process_active = true;
 
% Parameters not configurable by user 
    C = physconst('LightSpeed');
    % FMCW parameters
    active.PRF = 1/active.pulse_duration;
    active.slope = active.Bw/active.pulse_duration;
    active.F_Max = active.Bw/2;
    R_Max = beat2range(active.F_Max,active.slope);
    active.sample_duration = 1/active.Fs;
    active.active.samples_per_pulse = active.pulse_duration/active.sample_duration;
    active.number_pulses = capture_duration / active.pulse_duration;
    active.number_cap_samps = capture_duration/active.sample_duration;
    active.RF_freq = active.Fc/1e6;   % RF in MHz 
    active.Bw_M = active.Bw/1e6;      % BW in MHz
  

%% Create Sawtooth Chirp for bladeRF
chirp = saw_LFM_chirp(active.Bw,active.pulse_duration,active.Fs);
save_sc16q11('/tmp/chirp.sc16q11', chirp); %save chirp to binary file
clear chirp

% % %% Setup clock distribution
%  exit_code = setup_ref(1,3);
%  if exit_code > 0
%      "error configuring device"
%      return
%  end
%  pause(5);
%  exit_code = setup_ref(2,1);
%   if exit_code > 0
%      "error configuring device"
%      return
%  end
%  pause(5);
%  exit_code = setup_ref(3,1);
%  if exit_code > 0
%      "error configuring device"
%      return
%  end

%% Setup FMCW Radar
    % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''
    % Setup Tx SDR 
    [trig_flag_1,tx_command] = create_shell_command(true,...
                                   Experiment_ID,...
                                   active.number_cap_samps,... 
                                   active.number_pulses,...
                                   0,...
                                   active.Tx_gain,...
                                   active.Rx1_gain,...
                                   active.Rx2_gain,...
                                   active.RF_freq,...
                                   active.Bw_M,...
                                   active.Tx_SDR,...
                                   'slave',...
                                   3,...
                                   'tx');
    tx_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(tx_command);
    pause(5);


    % Setup Rx SDR 
    [trig_flag_2,rx_command] = create_shell_command(false,...
                                   Experiment_ID,...
                                   active.number_cap_samps,... 
                                   active.number_pulses,...
                                   0,...
                                   active.Tx_gain,...
                                   active.Rx1_gain,...
                                   active.Rx2_gain,...
                                   active.RF_freq,...
                                   active.Bw_M,...
                                   active.Rx_SDR,...
                                   'master',...
                                   1,...
                                   'rx'); 
    if trig_flag_1 && trig_flag_2
        "Trigger Conflict - FMCW Radar"
        return
    end
    %rx_command = rx_command + "&"; % uncomment for non-blocking system command execution                                                              
    system(rx_command); % Blocking system command execution    
    
    
 %% Save Raw Data and create header to directory 
    % make save directory
        make_dir = 'mkdir ' + exp_dir;
        system(make_dir);
    % move FMCW receive file to save directory
        move_file = 'mv /tmp/active_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
        rtn = system(move_file);
        if rtn == 0
            "FMCW Data Copied to Save directory"
        else 
            "FMCW Copy Failed"
            return

        end





%% FMCW Processing and Print RTI

    % load refsig for deramping
        refsig = load_refsig(active.Bw_M,active.Fs,active.Fc,active.pulse_duration);    
    % load Signal, Mix and Dermap Signal  
        file_location = exp_dir + 'active_' + Experiment_ID;
        lp_filter = true;
        [max_range_actual,deramped_signal,active.decimation_factor_actual] = deramp_and_decimate(file_location,active.max_range,refsig,capture_duration,active.number_pulses,active.Fs,active.slope,lp_filter);
        save(exp_dir + 'deramped_signal','deramped_signal')
    
    % Window and FFT Signal 
    % window signal
        w = window('blackman',size(deramped_signal,1));
        windowed_signal = deramped_signal.*w;
    % fft signal
        zero_padding = 1; % 1 = none; 2 = 100%
        processed_signal = fft(windowed_signal,size(windowed_signal,1)*zero_padding);
        beat_frequncies = processed_signal(1:(size(processed_signal,1)/2),:); % keep +ve beat frequencies

    % MTI Filtering 
        % Single Delay Line Filter 
            MTI_Data = zeros(size(beat_frequncies));
            active.range_bins = size(MTI_Data,1);
                  for i=2:active.number_pulses
                        MTI_Data(:,i) = beat_frequncies(:,i)-beat_frequncies(:,i-1);
                  end
            % IIR Filter
                [b, a] = butter(12, 0.04, 'high');
                  for i=1:active.range_bins
                        MTI_Data(i,:) = filtfilt(b,a,beat_frequncies(i,:));
                  end

    % Derive range and time axis 
        active.range_bins = 1:size(processed_signal,1);
        active.fftfrequncies =fftfreq(size(processed_signal,1),1/(active.Fs/active.decimation_factor_actual)); % possible beat frequencies
        active.slope = active.Bw/active.pulse_duration;
        ranges = (active.fftfrequncies*C)/(2*active.slope); % calculate true range bin size    
        active.range_axis = ranges(1:(size(ranges,2)/2)); % truncate to only +ve beat frequencies
        active.range_bin_size = ranges(2)
        active.time_axis = linspace(0,size(processed_signal,2)*active.pulse_duration,size(processed_signal,2));
      
     % Plot RTI
        RTI_plot= transpose(10*log10(abs(beat_frequncies./max(beat_frequncies(:)))));
        figure
        fig = imagesc(active.range_axis,active.time_axis,RTI_plot,[-50,0]);   
            ylabel('Time (s)')
            xlabel('Range (m)')
            title("FMCW RTI - " + Experiment_ID)
%             xlim([0 200])
            fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)    
    
     % Plot series of pulses
        fig = figure
        ranges_2_plot = floor(linspace(1,active.number_pulses,5));
        for i = ranges_2_plot
            plot(active.range_axis,RTI_plot(i,:));
            hold on
        end
            title("Single Pulse - " + Experiment_ID);
            xlim([0 200])
            grid on; grid minor;
            ylabel('Relative Power (dB)')
            xlabel('Range (m)')  
            fig_name = exp_dir + "Single_Pulse" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg') 

            
 % Spectrogram 
         % Parameters
            r_start = 5;
            r_stop = 12;
            l_fft = 256;
            pad_factor = 1;
            overlap_factor = 0.90;
     
         % Plot Spectrogram pre-MTI filtering
            integrated_data = sum(beat_frequncies(r_start:r_stop,:));
            [spect,active.doppler_axis] = spectrogram(integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,active.PRF,'centered','yaxis');
            spect= 10*log10(abs(spect./max(spect(:))));
            figure
            fig = imagesc(active.time_axis,-active.doppler_axis,spect,[-50 0]);   
                ylim([-500 500])
                c = colorbar
                c.Label.String='Norm Power (dB)'
                xlabel('Time (Sec)')
                % ylabel('Radial Velocity (mph)')   
                ylabel('Doppler Frequency (Hz)')  
                fig_title = "FMCW Spectrogram -" + Experiment_ID;
                title(fig_title);
                fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID + ".jpg";
                saveas(fig,fig_name,'jpeg')
                fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID;
                saveas(fig,fig_name)
                
          % Plot Spectrogram post-MTI filtering
            MTI_integrated_data = sum(MTI_Data(r_start:r_stop,:));
            [spect,f] = spectrogram(MTI_integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,active.PRF,'centered','yaxis');
            spect= 10*log10(abs(spect./max(spect(:))));
            figure
            fig = imagesc(active.time_axis,-f,spect,[-50 0]);   
                ylim([-500 500])
                c = colorbar
                c.Label.String='Norm Power (dB)'
                xlabel('Time (Sec)')
                % ylabel('Radial Velocity (mph)')   
                ylabel('Doppler Frequency (Hz)')  
                fig_title = "MTI FMCW Spectrogram - " + Experiment_ID;
                title(fig_title);
                fig_name = exp_dir + "MTI FMCW Spectrogram_" + Experiment_ID + ".jpg";
                saveas(fig,fig_name,'jpeg')
                fig_name = exp_dir + "MTI FMCW Spectrogram_" + Experiment_ID;
                saveas(fig,fig_name)

        %% Coherent integration 
            compressed_data = 10*log10(sum(beat_frequncies,2));
            figure
            plot(active.range_axis,abs(compressed_data))
            ylabel('Relative Power (dB)');xlabel('Range (m)');
            grid on; grid minor;

        %% Analyse phase of peak return 
            [value, pulse_rbin] = max(beat_frequncies(:,100)); 
            phase_9 = angle(processed_signal(pulse_rbin,:));
            phase_norm_9 = phase_9 - mean(phase_9(1:100));
            save(exp_dir + 'phase_norm_9', 'phase_9')
            figure
            plot(phase_norm_9)
            figure
            plot(abs(processed_signal(pulse_rbin,:)))
                    


