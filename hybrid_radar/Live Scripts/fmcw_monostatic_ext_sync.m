clear all
close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 122;    % Expeiment Name
capture_duration = 30;        % capture duration
save_directory = "~/Documents/bladeRAD_Captures/lab/"; % each experiment will save as a new folder in this directory
save_directory = "/media/sdrlaptop1/T7/08_08_2023_farm_n0/"; % each experiment will save as a new folder in this directory
exp_dir = save_directory + Experiment_ID + '/';


% FMCW Parameters 
active.Fs = 30e6;          % Sample Rate of SDR per I & Q (in reality Fs is double this)
active.pulse_duration = 0.5e-3;   % Desired Pulse Duration 
active.Bw = 30e6;          % LFM Bandwidth 
active.Fc = 2.44e9;   % Central RF 
active.Tx_gain = 66;       % [-23.75, 66] (S-Band = 23.5 dBm) (C-Band = 15.8 dBm)
active.Rx1_gain = 17;      % [-16, 60]
active.Rx2_gain = 0;       % [-16, 60]
Rx_1_lna = true;
Rx_2_lna = true;
active.Tx_SDR = 1;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
active.Rx_SDR = 2;   % SDR to use for RX
% Procesing Parameters
active.max_range = 1.4e5; %max range to LPF filter data to

process_active = false;
active_rd_slice = false;

 
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
% 
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
                                   2,...
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
                                   'slave',...
                                   1,...
                                   'rx'); 
    if trig_flag_1 && trig_flag_2
        "Trigger Conflict - FMCW Radar"
        return
    end
    %rx_command = rx_command + "&"; % uncomment for non-blocking system command execution                                                              
    system(rx_command); % Blocking system command execution 
    pause(5);
    

    
    
    
    
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
    save(exp_dir + 'Experimental Configuration') 





%% FMCW Processing and Print RTI
if process_active == true
    % load refsig for deramping
        refsig = load_refsig(active.Bw_M,active.Fs,active.Fc,active.pulse_duration);    
    % load Signal, Mix and Dermap Signal  
        file_location = exp_dir + 'active_' + Experiment_ID;
        lp_filter = true;
        [max_range_actual,deramped_signal,active.decimation_factor_actual] = deramp_and_decimate(file_location,active.max_range,refsig,capture_duration,active.number_pulses,active.Fs,active.slope,lp_filter);
        %save(exp_dir + 'deramped_signal','deramped_signal')
    
    % Window and FFT Signal 
    % window signal
        w = window('blackman',size(deramped_signal,1));
        deramped_signal = deramped_signal.*w;
    % fft signal
        zero_padding = 1; % 1 = none; 2 = 100%
        pre_processed_size = size(deramped_signal,1);
        deramped_signal = fft(deramped_signal,size(deramped_signal,1)*zero_padding);
        beat_frequncies = deramped_signal(1:(size(deramped_signal,1)/2),:); % keep +ve beat frequencies
        clear deramped_signal w

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
        active.range_bins = 1:pre_processed_size;
        active.fftfrequncies =fftfreq(pre_processed_size,1/(active.Fs/active.decimation_factor_actual)); % possible beat frequencies
        active.slope = active.Bw/active.pulse_duration;
        ranges = (active.fftfrequncies*C)/(2*active.slope); % calculate true range bin size    
        active.range_axis = ranges(1:(size(ranges,2)/2)); % truncate to only +ve beat frequencies
        active.range_bin_size = ranges(2);
        active.time_axis = linspace(0,size(beat_frequncies,2)*active.pulse_duration,size(beat_frequncies,2));
      
     % Plot RTI
        RTI_plot= transpose(20*log10(abs(beat_frequncies./max(beat_frequncies(:)))));
        figure;
        fig = imagesc(active.range_axis,active.time_axis,RTI_plot,[-50,0]);   
            ylabel('Time (s)');
            xlabel('Range (m)');
            title("FMCW RTI - " + Experiment_ID);
            xlim([0 400]);
            c = colorbar;
            c.Label.String='Norm Power [dB]';
            fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg');
            saveas(fig,fig_name);

    
     % Plot series of pulses
        fig = figure;
        ranges_2_plot = floor(linspace(1,active.number_pulses,5));
        for i = ranges_2_plot
            plot(active.range_axis,RTI_plot(i,:));
            xlim([0 400]);
            hold on
        end
            title("Single Pulse - " + Experiment_ID);
            xlim([0 400]);
            grid on; grid minor;
            ylabel('Relative Power (dB)');
            xlabel('Range (m)');  
            fig_name = exp_dir + "Single_Pulse" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg'); 
         clear RTI_plot
    % Cohernet Intgration
         coh_pulses = sum(transpose(beat_frequncies(:,1:10)));
         coh_pulses = coh_pulses./max(coh_pulses(:));
         singlepulse = abs(transpose(beat_frequncies(:,1)./max(beat_frequncies(:,1))));
         figure
         plot(20*log10(singlepulse)); hold on;
         plot(20*log10(abs(coh_pulses))); 
         grid on; grid minor; legend('Single Pulse','10 Integrated pulses')

            
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
            spect= 20*log10(abs(spect./max(spect(:))));
            figure
            fig = imagesc(active.time_axis,-active.doppler_axis,spect,[-50 0]);   
                %ylim([-500 500]);
                c = colorbar;
                c.Label.String='Norm Power (dB)';
                xlabel('Time (Sec)');
                % ylabel('Radial Velocity (mph)')   
                ylabel('Doppler Frequency (Hz)');  
                fig_title = "FMCW Spectrogram -" + Experiment_ID;
                title(fig_title);
                fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID + ".jpg";
                saveas(fig,fig_name,'jpeg')
                fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID;
                saveas(fig,fig_name)
                
          % Plot Spectrogram post-MTI filtering
            MTI_integrated_data = sum(MTI_Data(r_start:r_stop,:));
            [spect,f] = spectrogram(MTI_integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,active.PRF,'centered','yaxis');
            spect= 20*log10(abs(spect./max(spect(:))));
            figure
            fig = imagesc(active.time_axis,-f,spect,[-50 0]);   
                %ylim([-500 500])
                c = colorbar;
                c.Label.String='Norm Power (dB)';
                xlabel('Time (Sec)');
                % ylabel('Radial Velocity (mph)')   
                ylabel('Doppler Frequency (Hz)');  
                fig_title = "MTI FMCW Spectrogram - " + Experiment_ID;
                title(fig_title);
                fig_name = exp_dir + "MTI FMCW Spectrogram_" + Experiment_ID + ".jpg";
                saveas(fig,fig_name,'jpeg');
                fig_name = exp_dir + "MTI FMCW Spectrogram_" + Experiment_ID;
                saveas(fig,fig_name);

end

%% Process Active data into Range-Doppler Slices
if active_rd_slice == true 
           active.cpi = 0.5; % cohernet proccessing interval (s)
           active.cpi_overlap = 0.5; % overlap between CPIs (watch this - too large will cause slow exceution)
           active.zero_padding = 1; % 1 = none; 2 = 100%
           active.doppler_window = 'Hann';
           active.dynamic_range = 50;
           active.velocity_conv = C*(((1/C)/(active.Fc/C)));
                 active.range_bins = size(MTI_Data,1);
                 active.max_range_actual = max_range_actual; %1000;
  
          [active.number_cpi,....
           active.pulses_per_cpi,...
           active.range_doppler_slices] = rangeDopplerSlice(MTI_Data,active.cpi,active.PRF,...
                                                                 active.cpi_overlap,...
                                                                 active.zero_padding,...
                                                                 active.doppler_window);   
           
           % Determine Active Data Range and Doppler Axis                                              
             active.doppler_bins = 1:size(active.range_doppler_slices{1},1);
             active.doppler_axis = linspace(-active.PRF/2,active.PRF/2,size(active.doppler_bins,2));
             active.doppler_velocity_axis = active.doppler_axis * active.velocity_conv;
             active.range_axis = linspace(0,active.max_range_actual,active.range_bins);
             active.range_bins_axis = 1:active.range_bins;

           % Create video of range-Doppler slices
             video_name = exp_dir + "active_range-Doppler" + Experiment_ID + ".avi";
             %video_name = "range-Doppler_log_Exp_" + Experiment_ID + ".avi";       
             video_title = "Active Radar Capture";
             dynamic_range = +50;
             max_range = 100;
             max_doppler = +inf;
             frame_rate = 1/(capture_duration/active.number_cpi);    

             createVideo(active.range_doppler_slices,frame_rate,...
                         active.range_axis-10,max_range,...
                         -active.doppler_velocity_axis,max_doppler,...
                         dynamic_range,video_name,video_title);
end
