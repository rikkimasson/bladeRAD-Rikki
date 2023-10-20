 % paths to generic functions
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts/matlab/CFAR/',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/') 

% Local directory containing the radar data
    local_save_directory =  "/media/sdrlaptop1/T7/15_09_2022_farm/hybrid/";

% Processing Flags 
    process_active_flag = true;
    process_passive_flag = true;
    ground_truth_flag = true;
     
% Select the experiment you wish to process
    experiment_number = 2;

for i=experiment_number
    
%% Load .mat file containing experiment setup parameters
        mat_file_name = local_save_directory + i + "/Experimental Configuration.mat";
        load(mat_file_name);
        exp_dir = local_save_directory + i + '/';


%% Load Ground Truth Data 
if ground_truth_flag == true
    mat_file_name = local_save_directory  + i + "/geode1.mat";
    load(mat_file_name);

    gt.time_offset = 280-14; % offset between start time of radar and gps recording
    gt.range_offset = 8; % Range fudge factor between radar and gps
    gt.range = geode1.target_range(gt.time_offset:gt.time_offset+capture_duration) + gt.range_offset;
    gt.sog = geode1.SoG(gt.time_offset :gt.time_offset +capture_duration);

    % Interpolate to number of active radar pulses
        old_sample_points = 0:capture_duration;
        active_pulses_sample_points = 1:active.number_pulses;
        gt.pulses_interp_range = interp1(old_sample_points,gt.range,active_pulses_sample_points/1000,'linear');
end

%% FMCW Processing and Print RTI
if process_active_flag == true
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
                [b, a] = butter(10, 0.021, 'high');
                  for i=1:active.range_bins
                        MTI_Data(i,:) = filtfilt(b,a,beat_frequncies(i,:));
                  end

    % Derive range and time axis 
        active.n_range_bins = size(beat_frequncies,1);
        active.range_bins = 1:active.n_range_bins;
        active.fftfrequncies =fftfreq(size(processed_signal,1),1/(active.Fs/active.decimation_factor_actual)); % possible beat frequencies
        active.slope = active.Bw/active.pulse_duration;
        ranges = (active.fftfrequncies*C)/(2*active.slope); % calculate true range bin size    
        active.range_axis = ranges(1:(size(ranges,2)/2)); % truncate to only +ve beat frequencies
        active.range_bin_size = ranges(2)
        active.time_axis = linspace(0,size(processed_signal,2)*active.pulse_duration,size(processed_signal,2));
      
     % Plot RTI
        RTI_plot= transpose(20*log10(abs(beat_frequncies./max(beat_frequncies(:)))));
        figure
        fig = imagesc(active.range_axis,active.time_axis,RTI_plot,[-50,0]); hold on;
            if ground_truth_flag == true
                plot(gt.pulses_interp_range,active.time_axis,'-r', 'LineWidth', 2);
            end
            ylabel('Time (s)');xlabel('Range (m)');legend('Target Ground Truth')
            title("FMCW RTI - " + Experiment_ID)
            xlim([0 200])
            fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)

     % Plot MTI RTI
        RTI_plot= transpose(20*log10(abs(MTI_Data./max(MTI_Data(:)))));
        figure
        fig = imagesc(active.range_axis,active.time_axis,RTI_plot,[-50,0]); hold on;
           if ground_truth_flag == true
              plot(gt.pulses_interp_range,active.time_axis,'-r', 'LineWidth', 2);
           end
            ylabel('Time (s)');xlabel('Range (m)');legend('Target Ground Truth')
            title("FMCW RTI - " + Experiment_ID)
            xlim([0 200])
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

%% Process Active data into Range-Doppler Slices
           active.cpi = 0.5; % cohernet proccessing interval (s)
           active.cpi_overlap = 0.8; % overlap between CPIs (watch this - too large will cause slow exceution)
           active.zero_padding = 1; % 1 = none; 2 = 100%
           active.doppler_window = 'Hann';
           active.dynamic_range = 50
           active.velocity_conv = C*(((1/C)/(active.Fc/C)));
                 active.max_range_actual = max_range_actual; %1000;
          
           % Create range-Doppler surfaces  
              [active.number_cpi,....
               active.pulses_per_cpi,...
               active.range_doppler_slices] = rangeDopplerSlice(MTI_Data,active.cpi,active.PRF,...
                                                                     active.cpi_overlap,...
                                                                     active.zero_padding,...
                                                                     active.doppler_window);   
           
           % Create Active Data Range and Doppler Axis                                              
             active.doppler_bins = 1:size(active.range_doppler_slices{1},1);
             active.doppler_axis = linspace(-active.PRF/2,active.PRF/2,size(active.doppler_bins,2));
             active.doppler_velocity_axis = active.doppler_axis * active.velocity_conv;
             active.range_bins_axis = 1:active.range_bins;

           % Create video of range-Doppler slices
             video_name = exp_dir + "active_range-Doppler" + Experiment_ID + ".avi";
             video_name = "range-Doppler_log_Exp_" + Experiment_ID + ".avi";       
             video_title = "Active Radar Capture";
             dynamic_range = +100;
             max_range = 200;
             max_doppler = 60;
             frame_rate = 1/(capture_duration/active.number_cpi);    
             createVideo(active.range_doppler_slices,frame_rate,...
                         active.range_axis-10,max_range,...
                         -active.doppler_velocity_axis,max_doppler,...
                         dynamic_range,video_name,video_title);
                
end

 %% Passive Processing
 if process_passive_flag == true
    passive.max_range = 50; % maximum number of range bins - xcorr shifts.
    passive.range_zero_padding = 1; % 1 = none, 2 = 100%
    passive.td_corr = true; % true = time domain xcorr; false = freq domain xcorr
    passive.seg_s = 1000; % number of segments per second - analagos to PRF.
    passive.seg_percent = 100; % % of segment used for xcorr.
                            
    % load signal and split ref and sur
        file_location = exp_dir + 'passive_' + Experiment_ID;
        [ref_channel, sur_channel]  = load_passive_data(file_location);
    % Plot time domain signals
         figure
         fig = subplot(2,1,1);
            plot(real(ref_channel(1:4000000)));
            title("Ref channel time series");
         hold on
         subplot(2,1,2)
             plot(real(sur_channel(1:4000000)))
             title("Sur channel time series");    
             fig_name = exp_dir + "Time Domain Signals_" + Experiment_ID + ".jpg";
             saveas(fig,fig_name,'jpeg')

    % Add Noise to Passive Radar
        noisy_sur_channel= awgn(sur_channel,-2,'measured');

    % Batch process data and cross correlate
         [ref_matrix ,self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,noisy_sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.max_range,exp_dir,passive.range_zero_padding,passive.td_corr);
         save(exp_dir + 'passive_matrix','cc_matrix')

%     % RTI Plot
%         RTI_plot= transpose(10*log10(abs(cc_matrix./max(cc_matrix(:)))));
%         Range_bin = linspace(0,passive.max_range,size(cc_matrix,1));
%         time_axis = linspace(0,capture_duration,size(cc_matrix,2));
%         figure
%         fig = imagesc(Range_bin,time_axis,RTI_plot,[-50,0]);
%             grid on            
%             colorbar
%             ylabel('Time (Sec)')
%             xlabel('Range Bin')   
%             fig_title = "Passive RTI - " + Experiment_ID;
%             title(fig_title);
%             fig_name = exp_dir + "Passive RTI_" + Experiment_ID + ".jpg";
%             saveas(fig,fig_name,'jpeg'); saveas(fig,fig_name);
% 
%       % CAF of entire capture
%         f_axis = linspace(-passive.seg_s/2,passive.seg_s/2,size(cc_matrix,2));
%         t_cc_matrix = transpose(cc_matrix);
%         CAF = fftshift(fft(t_cc_matrix,size(t_cc_matrix,1),1),1);
%         figure
%         imagesc(Range_bin,f_axis,10*log10(abs(CAF./max(CAF(:)))),[-50 1]); 
%             ylim([-500 500])     
%             xlim([1 20])
%             c = colorbar
%             c.Label.String='Norm Power (dB)'
%             ylabel('Doppler Shift (Hz)')
%             xlabel('Range Bin')  
%             title("CAF for entire capture" + Experiment_ID)
%             fig_name = exp_dir + "CAF for entire capture_" + Experiment_ID + ".jpg";
%             saveas(fig,fig_name,'jpeg'); saveas(fig,fig_name);
%         
%      % Spectrogram 
%         int_bins = sum(cc_matrix(2:10,:),1);
%         r_bin = 1;
%         l_fft = 256;
%         pad_factor = 1;
%         overlap_factor = 0.99;
%         [spect,f] = spectrogram(int_bins,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,passive.seg_s,'centered','yaxis');
%         spect= 10*log10(abs(spect./max(spect(:))));
%         figure
%         fig = imagesc(time_axis,f,spect,[-30 0]);   
%             ylim([-500 +500])
%             c = colorbar
%             c.Label.String='Norm Power (dB)'
%             xlabel('Time (Sec)')
%             % ylabel('Radial Velocity (mph)')   
%             ylabel('Doppler Frequency (Hz)')  
%             fig_title = "Passive Spectrogram :- " + Experiment_ID;
%             title(fig_title);
%             fig_name = exp_dir + "Passive Spectrogram_" + Experiment_ID + ".jpg";
%             saveas(fig,fig_name,'jpeg'); saveas(fig,fig_name);
% 

%% Proccess Passive data into Range-Doppler Slices
           passive.PRF = passive.seg_s; % seg_s
           passive.cpi = 0.5; % coherent proccessing interval (s)
           passive.cpi_overlap = 0.8; % overlap between CPIs (watch this - too large will cause slow exceution)
           passive.doppler_window = 'Blackman-Harris';
           passive.dopp_zero_padding = 1;
           passive.dynamic_range = +50;
           passive.range_bin_size = (1/passive.Fs)/passive.range_zero_padding * C/2;   
           passive.max_range_m = passive.max_range*passive.range_bin_size;
           % Create range-Doppler surfaces 
                [passive.number_cpi,... 
                 passive.pulses_per_cpi,...
                 passive.range_doppler_slices] = rangeDopplerSlice(cc_matrix,passive.cpi,passive.PRF,...
                                                                   passive.cpi_overlap,...
                                                                   passive.dopp_zero_padding,...
                                                                   passive.doppler_window);  
           % Create Ambiguity surfaces  
                 [~,~,passive.self_ambg_slices] = rangeDopplerSlice(self_ambg_matrix,passive.cpi,passive.PRF,...
                                                                   passive.cpi_overlap,...
                                                                   passive.dopp_zero_padding,...
                                                                   passive.doppler_window);  
          
                          % Create Active Data Range and Doppler Axis                                              
  
            passive.cpi_stride = round(passive.pulses_per_cpi*(1-passive.cpi_overlap)); % number of pulses to stride each for next CPI
            passive.velocity_conv = C*(((1/C)/(passive.Fc/C)));
            passive.no_range_bins = size(cc_matrix,1);
            passive.range_bins = 1:passive.no_range_bins;
            passive.doppler_bins = passive.pulses_per_cpi*passive.dopp_zero_padding+1;
            passive.doppler_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.doppler_bins);
            passive.doppler_velocity_axis = passive.doppler_axis*passive.velocity_conv;
            passive.range_axis = linspace(0,passive.max_range_m,passive.no_range_bins);


            % create video of  range-Doppler slices
             video_name = exp_dir + "Range_Doppler_Slices" + Experiment_ID + ".avi";
             %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";       
             video_title = "Passive Pre-DSI";
             dynamic_range = +inf;
             max_range = 200;
             max_doppler = 100;
             frame_rate = 1/(capture_duration/passive.number_cpi);    
             createVideo(passive.range_doppler_slices,frame_rate,...
                         passive.range_axis,max_range,...
                         passive.doppler_axis,max_doppler,...
                         dynamic_range,video_name,video_title);
%                  
%% Direct Signal Interference Cancellation
   % set DSI cancellation parameters
      p = 0.999;        % subtraction parameter - P must be a positive 
                        % integer less than one to avoid unwanted discontinuities
                        % arising from zero values in the range–Doppler surface.
     threshold = 0.005; % cutoff threshold parameter
     max_iterations = 1000; % maximum number of itterations DSI is CLEANed from CAF slice
     number_rbins = size(passive.range_doppler_slices{1},2);
  
   % perform CLEAN based DSI Cancellation   
     passive.CLEANed_range_doppler_slices = CLEAN_DSI(passive.range_doppler_slices,...
                                               passive.self_ambg_slices,...
                                               ref_matrix,...
                                               passive.number_cpi,...
                                               passive.pulses_per_cpi,...
                                               passive.cpi_stride,...
                                               passive.dopp_zero_padding,...
                                               number_rbins,...
                                               max_iterations,...
                                               threshold,p,...
                                               passive.range_axis,passive.doppler_axis);

%   % create video of CLEANed range-Doppler slices
     video_name = exp_dir + "CLEANed_range-Doppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";
     %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";       
     video_title = "CLEANed Passive Radar Capture";
     dynamic_range = 50;
     max_range = 100;
     max_doppler = 20;
     frame_rate = 1/(capture_duration/passive.number_cpi);    
     createVideo(passive.CLEANed_range_doppler_slices,frame_rate,...
                 passive.range_axis,max_range,...
                 passive.doppler_velocity_axis,max_doppler,...
                 dynamic_range,video_name,video_title);

end


end
%% Interpolate ground truth to to number of CPI
       cpi_sample_points = linspace(0,capture_duration,active.number_cpi);
       gt.interp_range = interp1(old_sample_points,gt.range,cpi_sample_points,'linear');
       gt.interp_sog = interp1(old_sample_points,gt.sog,cpi_sample_points,'linear');
        
%% CFAR CAF Slices 

  %% Active CFAR
      % There are two CFAR detectors        
           % 1. Provides the detections cell array 'active.cfar.detections', 
           %    each cell is vector per CPI of detections in range and doppler.
           % 2. Provides cell array of binary valued range-Doppler surfaces 
           %    with detections as 1. 
     
      % Probability of false alarm
        pfa = 1e-12;
      
      % 1. Array of vectors with detection in range and Doppler
            % Create the CFAR detector.
                active.detector = mod_CFARDetector2D('TrainingBandSize',[4,4], ...
                                                 'GuardBandSize',[3,3],...
                                                 'ThresholdFactor','Auto',...
                                                 'Method','CA', ...
                                                 'ProbabilityFalseAlarm',pfa,...
                                                 'ThresholdOutputPort',true,...
                                                 'NoisePowerOutputPort',true,...
                                                 'OutputFormat','Detection index');             
             
            % Create CUT Matrix
                active.full_cutidx = createCutMatrix(active.detector,active.n_range_bins,size(active.doppler_axis,2));                            
      
            % loop through matrix and operate CFAR to get detections, threshold and noise power estimate    
                [active.cfar.detections,...
                 active.cfar.threshold,...
                 active.cfar.noise_estimate] = cfarSlicesDetections(active.range_doppler_slices,active.detector,...
                                                    active.full_cutidx,active.range_axis,active.doppler_velocity_axis);
            % plot detections in range and Doppler
                figure
                for i=1:active.number_cpi
                    plot(active.cfar.detections{i}(1,:),active.cfar.detections{i}(2,:),'r*')
                    hold on
                end
                grid on; grid minor; ylabel('Doppler Shift (Hz)'); xlabel('Range (m)'); xlim([-inf 200]) 

%       % 2. Array of vectors with detection in range and Doppler
%         % Create the CFAR detector.
%             active.detector = mod_CFARDetector2D('TrainingBandSize',[4,4], ...
%                                              'GuardBandSize',[3,3],...
%                                              'ThresholdFactor','Auto',...
%                                              'Method','CA', ...
%                                              'ProbabilityFalseAlarm',pfa,...
%                                              'ThresholdOutputPort',true);
% 
%          % loop through matrix and operate CFAR to get range-Doppler CFAR results   
%              [active.cfar.range_doppler_slices] = cfarSlices(active.range_doppler_slices,active.detector,active.full_cutidx);
%   
%          % create video of CFAR data
%              video_name = exp_dir + "CFAR_active_range-Doppler_Exp_" + Experiment_ID + ".avi";
%              video_title = "CFAR Active Radar Capture";
%              max_range = +150;
%              max_doppler = +inf;
%              frame_rate = 1/(capture_duration/active.number_cpi);   
%              createCFARVideo(active.cfar.range_doppler_slices,frame_rate,...
%                              active.range_axis,max_range,...
%                             -active.doppler_axis,max_doppler,...
%                             video_name,video_title);
    
    %% CFAR Passive Data 
     % There are two CFAR detectors        
           % 1. Provides the detections cell array 'active.cfar.detections', 
           %    each cell is vector per CPI of detections in range and doppler.
           % 2. Provides cell array of binary valued range-Doppler surfaces 
           %    with detections as 1.

     % 1. Array of vectors with detection in range and Doppler   
        % Create the CFAR detector.
            passive.detector = mod_CFARDetector2D('TrainingBandSize',[4,4], ...
                                             'GuardBandSize',[5,5],...
                                             'ThresholdFactor','Auto',...
                                             'Method','CA', ...
                                             'ProbabilityFalseAlarm',pfa,...
                                             'ThresholdOutputPort',true,...
                                             'NoisePowerOutputPort',true,...
                                             'OutputFormat','Detection index');             
       % Create CUT Matrix
            passive.full_cutidx = createCutMatrix(passive.detector,passive.no_range_bins,size(passive.doppler_axis,2));                            
 
       % Loop through matrix and operate CFAR to get detections, threshold and noise power estimate    
            [passive.cfar.detections,...
             passive.cfar.threshold,...
             passive.cfar.noise_estimate] = cfarSlicesDetections(passive.CLEANed_range_doppler_slices,passive.detector,...
                                                passive.full_cutidx,passive.range_axis,-passive.doppler_velocity_axis);
         
       % Shift passive radar range to align with active radar range
             passive.range_shift = 2 * active.range_bin_size;
            
       % Plot active and passive detections in range and Doppler
            figure
            plot(gt.interp_range,2*-gt.interp_sog,'black-',LineWidth=2);hold on;      
            for i=1:active.number_cpi
                plot(active.cfar.detections{i}(1,:),active.cfar.detections{i}(2,:),'r*');hold on;        
                plot(passive.cfar.detections{i}(1,:)+passive.range_shift,passive.cfar.detections{i}(2,:),'b*');hold on;    
            end
            grid on; ylabel('Doppler Velocity (m/s)'); xlabel('Range (m)'); legend('Active Radar Detections','Passive Radar Detections');
            xlim([-inf 100]); ylim([-20 20]);
            
%       % 2. Array of vectors with detection in range and Doppler
%         % Create the CFAR detector.
%             passive.detector = mod_CFARDetector2D('TrainingBandSize',[4,4], ...
%                                              'GuardBandSize',[3,3],...
%                                              'ThresholdFactor','Auto',...
%                                              'Method','CA', ...
%                                              'ProbabilityFalseAlarm',pfa,...
%                                              'ThresholdOutputPort',true);
%         
%         % loop through matrix and operate CFAR to get range-Doppler CFAR results   
%              [passive.cfar.range_doppler_slices] = cfarSlices(passive.CLEANed_range_doppler_slices,passive.detector,passive.full_cutidx);
%        
%         % create video of CFAR data
%              video_name = exp_dir + "CFAR_passvie_range-Doppler_Exp_" + Experiment_ID + ".avi";
%              video_title = "CFAR Active Radar Capture";
%              max_range = 200;
%              max_doppler = 200;
%              frame_rate = 1/(capture_duration/passive.number_cpi);    
%              figure;
%              createCFARVideo(passive.cfar.range_doppler_slices,frame_rate,...
%                              passive.range_axis,max_range,...
%                             -passive.doppler_axis,max_doppler,...
%                             video_name,video_title);

%% Hybrid Radar Signal Processing

    %% Interpolate passive range-Doppler surfaces to same cell size as active radar 
        % Interplolation factors 
            range_factor = passive.range_bin_size/active.range_bin_size; 
            doppler_factor = passive.velocity_conv/active.velocity_conv;  
        % Interpolate passive range-Doppler surfaces
            passive.interp_range_doppler_slices = interpolateRangeDoppler(passive.CLEANed_range_doppler_slices,range_factor,doppler_factor);
        % Interpolated axis
            [passive.interp_doppler_bins, passive.no_interp_range_bins] = size(passive.interp_range_doppler_slices{1});
            passive.interp_doppler_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.interp_doppler_bins);
            passive.interp_doppler_velocity_axis = passive.interp_doppler_axis*passive.velocity_conv;
            passive.interp_range_axis = linspace(0,passive.max_range_m,passive.no_interp_range_bins);
            passive.interp_range_bins = 1:passive.interp_no_range_bins;

% %     %% Create Hybrid Video of both channels of data 
%        video_name = exp_dir + "hybrid_range-Doppler_Exp_" + Experiment_ID + ".avi";     
%        dynamic_range = 50;
%        range_limit = 200;
%        doppler_limit = 20;
%        createVideos(active.range_doppler_slices,...
%                      passive.interp_range_doppler_slices,...
%                      frame_rate,...
%                      active.range_axis,...
%                      passive.interp_range_axis,...
%                      -active.doppler_velocity_axis,...
%                      passive.interp_doppler_velocity_axis,...
%                      range_limit,...
%                      doppler_limit,...
%                      dynamic_range,...
%                      video_name,video_title);
                                        

    %% CFAR Interpolated Passive Data
    pfa = 1e-12;     
    % There are two CFAR detectors        
               % 1. Provides the detections cell array 'active.cfar.detections', 
               %    each cell is vector per CPI of detections in range and doppler.
               % 2. Provides cell array of binary valued range-Doppler surfaces 
               %    with detections as 1.
    
         % 1. Array of vectors with detection in range and Doppler   
            % Create the CFAR detector.
                passive.detector = mod_CFARDetector2D('TrainingBandSize',[8,8], ...
                                                 'GuardBandSize',[10,10],...
                                                 'ThresholdFactor','Auto',...
                                                 'Method','CA', ...
                                                 'ProbabilityFalseAlarm',pfa,...
                                                 'ThresholdOutputPort',true,...
                                                 'NoisePowerOutputPort',true,...
                                                 'OutputFormat','Detection index');             
           % Create CUT Matrix
                passive.interp_full_cutidx = createCutMatrix(passive.detector,passive.interp_no_range_bins,passive.interp_doppler_bins);                            
     
           % Loop through matrix and operate CFAR to get detections, threshold and noise power estimate    
                [passive.cfar.interp_detections,...
                 passive.cfar.interp_threshold,...
                 passive.cfar.interp_noise_estimate] = cfarSlicesDetections(passive.interp_range_doppler_slices,passive.detector,...
                                                     passive.interp_full_cutidx,passive.interp_range_axis,-passive.interp_doppler_velocity_axis);
             
           % Shift passive radar range to align with active radar range
                 passive.range_shift = 2 * active.range_bin_size;
                
           % Plot active and passive detections in range and Doppler
                figure
                plot(gt.interp_range-10,2*-gt.interp_sog,'black-',LineWidth=2);hold on;      
                for i=1:active.number_cpi
                    plot(active.cfar.detections{i}(1,:)-10,active.cfar.detections{i}(2,:),'r*');hold on;        
                    plot(passive.cfar.interp_detections{i}(1,:)-10+passive.range_shift,passive.cfar.interp_detections{i}(2,:),'b*');hold on;    
                end
                grid on; ylabel('Doppler Velocity (m/s)'); xlabel('Range (m)'); legend('Ground Truth','Active Radar Detections','Passive Radar Detections');
                xlim([0 100]); ylim([-20 20]);
                
          % 2. Array of vectors with detection in range and Doppler
            % Create the CFAR detector.
                passive.detector = mod_CFARDetector2D('TrainingBandSize',[4,4], ...
                                                 'GuardBandSize',[3,3],...
                                                 'ThresholdFactor','Auto',...
                                                 'Method','CA', ...
                                                 'ProbabilityFalseAlarm',pfa,...
                                                 'ThresholdOutputPort',true);
            
            % loop through matrix and operate CFAR to get range-Doppler CFAR results   
                 [passive.cfar.interp_range_doppler_slices] = cfarSlices(passive.interp_range_doppler_slices,passive.detector,passive.full_cutidx);
           
       
        % create video of CFAR data
             video_name = exp_dir + "CFAR_passvie_range-Doppler_Exp_" + Experiment_ID + ".avi";
             video_title = "CFAR Active Radar Capture";
             max_range = 200;
             max_doppler = 200;
             frame_rate = 1/(capture_duration/passive.number_cpi);    
             figure;
             createCFARVideo(passive.cfar.interp_range_doppler_slices,frame_rate,...
                             passive.range_axis,max_range,...
                            -passive.doppler_axis,max_doppler,...
                            video_name,video_title);

       
%% Plot range and Doppler of peak cell in active and passive range-Doppler surfaces.
  % Create empty vectors for detections 
    passive_range_detections = zeros(active.number_cpi,1);
    active_range_detections = zeros(active.number_cpi,1);
    passive_doppler_detections = zeros(active.number_cpi,1);
    active_doppler_detections = zeros(active.number_cpi,1);
  % Create lookup matrix to find range and doppler of detection on matrix     
    passive_range_lookup = ones(size(passive.range_doppler_slices{1})).*passive.range_axis;
    passive_doppler_lookup = ones(size(passive.range_doppler_slices{1})).*transpose(passive.doppler_velocity_axis);
    passive_range_lookup = ones(size(passive.interp_range_doppler_slices{1})).*passive.interp_range_axis;
    passive_doppler_lookup = ones(size(passive.interp_range_doppler_slices{1})).*transpose(passive.interp_doppler_velocity_axis);
    active_range_lookup = ones(size(active.range_doppler_slices{1})).*active.range_axis;
    active_doppler_lookup = ones(size(active.range_doppler_slices{1})).*transpose(active.doppler_velocity_axis);
  % Find peak cell on range doppler surface for each CPI
    for j=1:active.number_cpi   
       
%        [mag, index] = max(passive.CLEANed_range_doppler_slices{j}(:));
%        passive_range_detections(j) = passive_range_lookup(index);
%        passive_doppler_detections(j) = passive_doppler_lookup(index);
        
       [mag, index] = max(passive.interp_range_doppler_slices{j}(:));
       passive_range_detections(j) = passive_range_lookup(index);
       passive_doppler_detections(j) = passive_doppler_lookup(index);
    
       [mag, index] = max(active.range_doppler_slices{j}(:));
       active_range_detections(j) = active_range_lookup(index);
       active_doppler_detections(j) = active_doppler_lookup(index);         
    end

    passive.range_shift = 2 * passive.interp_range_axis(2);
    % Plot peak detection in range over time   
        time_axis = linspace(0,capture_duration,active.number_cpi);
        fig = figure
        pa = plot(time_axis,passive_range_detections+passive.range_shift); hold on;
        a = plot(time_axis,active_range_detections);hold on;
        c = plot(time_axis,abs(active_range_detections-(passive_range_detections+passive.range_shift)),'-g');
        pa.Marker = '*'; a.Marker = '*';
        plot(0:capture_duration,gt.range,'black-')
        grid on; grid minor; ylabel('Range (m)'); xlabel('Time (s)');
        % ylim([-inf 50])
        legend('Passive Radar Peak Return','Active Radar Peak Return','Relative Sensor Error','GPS Ground Truth')
        fig_name = exp_dir + "Hybrid Data RT - " + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg'); saveas(fig,fig_name);
    
        std((active_range_detections(7:end)-(passive_range_detections(7:end)+passive.range_shift)))

    % Plot peak detection in velocity over time   
        fig = figure
        plot(time_axis,-passive_doppler_detections,'-bl*'); hold on;
        a = plot(time_axis,active_doppler_detections,'-r*'); hold on;
        c = plot(time_axis,active_doppler_detections+passive_doppler_detections,'g')
        plot(0:capture_duration,2*-gt.sog,'black-')
        grid on; grid minor;
        % ylim([-8 8])
        ylabel('Velocity (m/s)'); xlabel('Time (s)');
        legend('Passive Radar Peak Return','Active Radar Peak Return','Relative Sensor Error','GPS Ground Truth')
        fig_name = exp_dir + "Hybrid Data DT - " + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg'); saveas(fig,fig_name);
    
        mean((active_doppler_detections(7:end)+(passive_doppler_detections(7:end))))

%% Calculate SNR of peak return 

    passive_detection_snr = zeros(active.number_cpi,1);
    active_detection_snr = zeros(active.number_cpi,1);
   % loop thorugh CPIs 
    for j=1:active.number_cpi   
       [target_mag, index] = max(passive.interp_range_doppler_slices{j}(:));
       passive_range_detections(j) = passive_range_lookup(index)
       passive_target_power = 10*log10(abs(target_mag))
       passive_noise_estimate = 10*log10(mean(abs(passive.interp_range_doppler_slices{j}(:,110:130)),'all'))   
       passive_detection_snr(j) = passive_target_power - passive_noise_estimate;
     
       [target_mag, index] = max(active.range_doppler_slices{j}(:));
       active_range_detections(j) = active_range_lookup(index);
       active_target_power = 10*log10(abs(target_mag));
       active_noise_estimate = 10*log10(mean(abs(active.range_doppler_slices{j}(:,80:100)),'all'));
       active_detection_snr(j) = active_target_power - active_noise_estimate;
    end  
    
    fig = figure    
    plot(passive_range_detections,passive_detection_snr,'*r'); hold on;
    plot(active_range_detections,active_detection_snr,'*b')    
    legend('SNR Passive Radar Peak Return','SNR Active Radar Peak Return')
    grid on; grid minor;
    ylabel('SNR (dB)'); xlabel('Range (m)');
    % xlim([-inf 100])
    fig_name = exp_dir + "Passive SNR of Peak Detection - CPI = " + passive.cpi + " s.jpg";
    saveas(fig,fig_name,'jpeg'); saveas(fig,fig_name);


%% Plot Active and Passive CFAR detections + Ground Truth 

% plot detections and ground truth in range and Doppler
        figure
        plot(gt.interp_range,2*-gt.interp_sog,'black-','LineWidth',2)
        hold on
        for i=1:active.number_cpi
            plot(active.cfar.detections{i}(1,:)-2,active.cfar.detections{i}(2,:),'r*')
            hold on
            plot(passive.cfar.detections{i}(1,:)+8,passive.cfar.detections{i}(2,:),'b*')
            hold on
        end
            grid on; grid minor
            ylabel('Doppler Velocity (m/s)')
            xlabel('Range (m)')  
            xlim([-inf 100])
            ylim([-20 20])
            legend('GPS Ground Truth','Active Radar Detections','Passive Radar Detections')

%% Calculate Scenario Pd
    % Currently 'computePD' function is setup just to check if there is a 
    % detection at the correct range, not range and Doppler.
    
    % compute the active sensor probability of detection
        active.cfar.range_error_lim = 2.5; 
        active.cfar.sog_error_lim = 200;
        active.rangeBias = -2;
        [active.cfar.detected, active.cfar.pd] = computePD(active.cfar.detections,gt.interp_range,2*-gt.interp_sog,active.cfar.range_error_lim,active.cfar.sog_error_lim,active.rangeBias);
        figure; plot(active.cfar.detected); ylim([0 2]);title('Active Detections');active.cfar.pd
        "Active Pd : " + (active.cfar.pd)

    % compute the passive sensor probability of detection
        passive.cfar.range_error_lim = 2.5; 
        passive.cfar.sog_error_lim = 200;
        passive.rangeBias = passive.range_shift;
        [passive.cfar.detected, passive.cfar.pd] = computePD(passive.cfar.interp_detections,gt.interp_range,2*-gt.interp_sog,passive.cfar.range_error_lim,passive.cfar.sog_error_lim,passive.rangeBias);
        figure; plot(passive.cfar.detected); ylim([0 2]); title('Passive Detections'); 
        "Passive Pd : " + (passive.cfar.pd)
         

    % fused detection result
        fused = active.cfar.detected + passive.cfar.detected;
        fused = ceil(fused/2);
        fused_pd = (sum(fused)/active.number_cpi)*100;
        figure; plot(fused); ylim([0 2]);title('Fused Detections');
        "Hybrid Pd : " + (fused_pd)
    

