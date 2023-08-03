clear all
% close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 1005;    % Expeiment Name
capture_duration = 30;        % capture duration
save_directory = "~/Documents/bladeRAD_Captures/lab/"; % each experiment will save as a new folder in this directory
exp_dir = save_directory + Experiment_ID + '/';


% FMCW Parameters 
active.Fs = 30e6;          % Sample Rate of SDR per I & Q (in reality Fs is double this)
active.pulse_duration = 0.2e-3;   % Desired Pulse Duration 
active.Bw = 30e6;          % LFM Bandwidth 
active.Fc = 2.44e9;   % Central RF 
active.Tx_gain = 66;       % [-23.75, 66] (S-Band = 23.5 dBm) (C-Band = 15.8 dBm)
active.Rx1_gain = 3;      % [-16, 60]
active.Rx2_gain = 0;       % [-16, 60]
Rx_1_lna = true;
Rx_2_lna = true;
active.Tx_SDR = 1;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
active.Rx_SDR = 2;   % SDR to use for RX
% Procesing Parameters
active.max_range = 1.4e5; %max range to LPF filter data to

% Passive Radar Parameters
passive.Fc = 690e6;   % Central RF    sample_duration
passive.Ref_gain = 0; % 26dB seems good for C-Band Patch Antennas [-16, 60]
passive.Sur_gain = 26; % 40dB seems good for C-Band Patch Antennas [-16, 60]
passive.SDR = 3;   % SDR to use for Passive Radar - labelled on RFIC Cover and bladeRAD Facia Panel
passive.Bw = 10e6;
passive.Fs = passive.Bw;
passive.max_range = 1000; %max range to cross-correlate to

process_active = false;
active_rd_slice = false;
process_passive = false;
 
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
    
    % passive parameters
    passive_Fs = passive.Bw;
    passive_sample_duration = 1/passive_Fs;
    passive_number_cap_samps = (2*capture_duration)/passive_sample_duration;
    passive_RF_freq = passive.Fc/1e6;   % RF in MHz 
    passive.Bw_M = passive.Bw/1e6;      % BW in MHz


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
    [trig_flag_2,rx_command] = create_shell_command(true,...
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
    rx_command = rx_command + "&"; % uncomment for non-blocking system command execution                                                              
    system(rx_command); % Blocking system command execution 
    pause(5);
    
% %% Setup Passive Radar
    % by default clock out is always enabled 
    % 1 'set clock_sel external'; 
    % 2 'set clock_ref enable';

    % Setup Passive SDR 
   [trig_flag_3, passive_command] = create_shell_command(false,...
                                   Experiment_ID,...
                                   passive_number_cap_samps,... 
                                   0,...
                                   0,...
                                   0,...
                                   passive.Ref_gain,...
                                   passive.Sur_gain,...
                                   passive_RF_freq,...
                                   passive.Bw_M,...
                                   passive.SDR,...
                                   'slave',...
                                   1,...
                                   'pass');
                               
    if trig_flag_1 && trig_flag_3 || trig_flag_2 && trig_flag_3
        "Trigger Conflict - Passive"
        return
    end
%     passive_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(passive_command);
    pause(2)





    
    
    
    
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
    % move passive file to save directory
    move_file = 'mv /tmp/passive_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Passive Data Copied to Save directory"
    else 
        "Passive Copy Failed"
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
            xlim([0 200]);
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
            xlim([0 200]);
            hold on
        end
            title("Single Pulse - " + Experiment_ID);
            xlim([0 200]);
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

 %% Passive Processing
 if process_passive == true
    passive.max_range = 50
    passive.range_zero_padding = 1; % 1 = none, 2 = 100%
    passive.td_corr = true; % if True Cross correlation done in the time domain; if flase FD correlation
    passive.seg_s = 1000; % number of segments per second - analagos to PRF.
    passive.seg_percent = 100;  % percentage of segment used for cross coreclation of 
                            % survallance and reference. Will affect SNR dramatically.
    % load signal and split ref and sur
        file_location = exp_dir + 'passive_' + Experiment_ID;
        [ref_channel, sur_channel]  = load_passive_data(file_location);
    % Plot time domain signals
         figure
         fig = subplot(2,1,1);
            plot(real(ref_channel(1:4000000)));
            title("Ref channel time series: " + max(abs(real(ref_channel))));
         hold on
         subplot(2,1,2)
             plot(real(sur_channel(1:4000000)))
             title("Sur channel time series: " + max(abs(real(sur_channel))));    
             fig_name = exp_dir + "Time Domain Signals_" + Experiment_ID + ".jpg";
             saveas(fig,fig_name,'jpeg')
             
    % Batch process data and cross correlate
         [ref_matrix ,self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.max_range,exp_dir,passive.range_zero_padding,passive.td_corr);
         save(exp_dir + 'passive_matrix','cc_matrix')
         clear ref_channel sur_channel

    % RTI Plot
        RTI_plot= transpose(20*log10(abs(cc_matrix./max(cc_matrix(:)))));
        Range_bin = linspace(0,passive.max_range,size(cc_matrix,1));
        time_axis = linspace(0,capture_duration,size(cc_matrix,2));
        figure
        fig = imagesc(Range_bin,time_axis,RTI_plot,[-50,0]); 
            %xlim([1 20])
            grid on            
            colorbar
            ylabel('Time (Sec)')
            xlabel('Range Bin')   
            fig_title = "Passive RTI - " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "Passive RTI_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)
         clear RTI_plot

      % CAF of entire capture
        f_axis = linspace(-passive.seg_s/2,passive.seg_s/2,size(cc_matrix,2));
        t_cc_matrix = transpose(cc_matrix);
        CAF = fftshift(fft(t_cc_matrix,size(t_cc_matrix,1),1),1);
        figure
        imagesc(Range_bin,f_axis,20*log10(abs(CAF./max(CAF(:)))),[-50 1]); 
            ylim([-500 500])     
            xlim([1 20])
            c = colorbar
            c.Label.String='Norm Power (dB)'
            ylabel('Doppler Shift (Hz)')
            xlabel('Range Bin')  
            title("CAF for entire capture" + Experiment_ID)
            fig_name = exp_dir + "CAF for entire capture_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)
        clear CAF
        

     % Spectrogram 
        int_bins = sum(cc_matrix(2:10,:),1);
        r_bin = 1;
        l_fft = 256;
        pad_factor = 1;
        overlap_factor = 0.99;
        [spect,f] = spectrogram(int_bins,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,passive.seg_s,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1Blackman-Harris:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/passive.Fc)*2.237;
        spect= 20*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-30 0]);   
            ylim([-500 +500])
            c = colorbar
            c.Label.String='Norm Power (dB)'
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "Passive Spectrogram :- " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "Passive Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)

%% Proccess Passive data into Range-Doppler Slices
           passive.PRF = passive.seg_s; %seg_s
           passive.cpi = 0.5; % cohernet proccessing interval (s)
           passive.cpi_overlap = 0.5; % overlap between CPIs (watch this - too large will cause slow exceution)
           passive.doppler_window = 'Hann';
           passive.zero_padding = 1;
           passive.dynamic_range = +inf;
           passive.max_range = passive.max_range;       
     
                [passive.number_cpi,... 
                 passive.pulses_per_cpi,...
                 passive.range_doppler_slices] = rangeDopplerSlice(cc_matrix,passive.cpi,passive.PRF,...
                                                                   passive.cpi_overlap,...
                                                                   passive.zero_padding,...
                                                                   passive.doppler_window);  
                
                 [~,~,passive.self_ambg_slices] = rangeDopplerSlice(self_ambg_matrix,passive.cpi,passive.PRF,...
                                                                   passive.cpi_overlap,...
                                                                   passive.zero_padding,...
                                                                   passive.doppler_window);  
            
            passive.cpi_stride = round(passive.pulses_per_cpi*(1-passive.cpi_overlap)); % number of pulses to stride each for next CPI
            passive.velocity_conv = C*(((1/C)/(passive.Fc/C)));
            passive.range_bins = size(cc_matrix,1);
            passive.doppler_bins = passive.pulses_per_cpi*passive.zero_padding+1;
            passive.doppler_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.doppler_bins);
            passive.doppler_velocity_axis = passive.doppler_axis*passive.velocity_conv;
            passive.range_axis = linspace(0,passive.max_range,passive.range_bins);


            % create video of  range-Doppler slices
             video_name = exp_dir + "Range_Doppler_Slices" + Experiment_ID + ".avi";
             %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";       
             video_title = "Passive Pre-DSI";
             dynamic_range = +inf;
             max_range = 10;
             max_doppler = 100;
             frame_rate = 1/(capture_duration/passive.number_cpi);    
             createVideo(passive.range_doppler_slices,frame_rate,...
                         passive.range_axis,max_range,...
                         passive.doppler_axis,max_doppler,...
                         dynamic_range,video_name,video_title);
        
                 
%% Direct Signal Interference Cancellation
   % set DSI cancellation parameters
      p = 0.999;        % subtraction parameter - P must be a positive 
                        % integer less than one to avoid unwanted discontinuities
                        % arising from zero values in the range–Doppler surface.
     threshold = 0.005; % cutoff threshold parameter
     max_iterations = 100; % maximum number of itterations DSI is CLEANed from CAF slice
     number_rbins = size(passive.range_doppler_slices{1},2);
  
   % perform CLEAN based DSI Cancellation   
     passive.CLEANed_range_doppler_slices = CLEAN_DSI(passive.range_doppler_slices,...
                                               passive.self_ambg_slices,...
                                               ref_matrix,...
                                               passive.number_cpi,...
                                               passive.pulses_per_cpi,...
                                               passive.cpi_stride,...
                                               passive.zero_padding,...
                                               number_rbins,...
                                               max_iterations,...
                                               threshold,p,...
                                               passive.range_axis,passive.doppler_axis);

close all
  % create video of CLEANed range-Doppler slices
     video_name = exp_dir + "CLEANed_range-Doppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";
     %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";       
     video_title = "CLEANing Passive Radar Capture";
     dynamic_range = 50;
     max_range = 20;
     max_doppler = 50;
     frame_rate = 1/(capture_duration/passive.number_cpi);    
     createVideo(passive.CLEANed_range_doppler_slices,frame_rate,...
                 passive.range_axis,max_range,...
                 passive.doppler_axis,max_doppler,...
                 dynamic_range,video_name,video_title);

 end