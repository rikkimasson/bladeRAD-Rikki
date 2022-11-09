% clear all
% % close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts/matlab/CFAR/',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions
        
save_directory1 =  "/media/piers/T7/15_09_2022_farm/hybrid/";

process_active_a = true;
process_passive_a = true;
 
%% create array with experiment names
dinfo = dir(save_directory1);
names_cell = {dinfo.name};
names_cell = names_cell(3:end);

for i=2
    
%% load .mat file containing experiment parameters
        mat_file_name = save_directory1 + i + "/Experimental Configuration.mat";
        load(mat_file_name);
   
        exp_dir = save_directory1 + i + '/';

%% FMCW Processing and Print RTI
if process_active_a == true
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
        fig = imagesc(active.range_axis,active.time_axis,RTI_plot,[-50,0]); hold on;
        plot(interp_flight_path,active.time_axis,'-r', 'LineWidth', 2);
            ylabel('Time (s)')
            xlabel('Range (m)')
            title("FMCW RTI - " + Experiment_ID)
            %xlim([0 200])
            fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)    
    


    grid on            
    c = colorbar
    c.Label.String='Norm Power (dB)'
    ylabel('Time (Sec)')
    xlabel('Range (m)')            
    xlim([0 200])
    legend('Target GPS Ground Truth')








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

     
      %% Spectrogram 
        r_start = 5;
        r_stop = 12;
        l_fft = 256;
        pad_factor = 1;
        overlap_factor = 0.90;
        
        integrated_data = sum(processed_signal(r_start:r_stop,:));
        [spect,active.doppler_axis] = spectrogram(integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,active.PRF,'centered','yaxis');
%         v=dop2speed(f,C/passive.Fc)*2.237;
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
            
            MTI_integrated_data = sum(MTI_Data(r_start:r_stop,:));
            [spect,f] = spectrogram(MTI_integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,active.PRF,'centered','yaxis');
            v=dop2speed(f,C/passive.Fc)*2.237;
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
             max_doppler = 20;
             frame_rate = 1/(capture_duration/active.number_cpi);    

             createVideo(active.range_doppler_slices,frame_rate,...
                         active.range_axis-10,max_range,...
                         -active.doppler_velocity_axis,max_doppler,...
                         dynamic_range,video_name,video_title);




                
end

 %% Passive Processing
 if process_passive_a == true
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
            title("Ref channel time series");
         hold on
         subplot(2,1,2)
             plot(real(sur_channel(1:4000000)))
             title("Sur channel time series");    
             fig_name = exp_dir + "Time Domain Signals_" + Experiment_ID + ".jpg";
             saveas(fig,fig_name,'jpeg')
    
 
             
    % Batch process data and cross correlate
         [ref_matrix ,self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.max_range,exp_dir,passive.range_zero_padding,passive.td_corr);
         save(exp_dir + 'passive_matrix','cc_matrix')

    % RTI Plot
        RTI_plot= transpose(10*log10(abs(cc_matrix./max(cc_matrix(:)))));
        Range_bin = linspace(0,passive.max_range,size(cc_matrix,1));
        time_axis = linspace(0,capture_duration,size(cc_matrix,2));
        figure
        fig = imagesc(Range_bin,time_axis,RTI_plot,[-50,0]);
%             xlim([1 20])
            grid on            
            colorbar
            ylabel('Time (Sec)')
            xlabel('Range Bin')   
            fig_title = "Passive RTI - " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "Passive RTI_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)


      % CAF of entire capture
        f_axis = linspace(-passive.seg_s/2,passive.seg_s/2,size(cc_matrix,2));
        t_cc_matrix = transpose(cc_matrix);
        CAF = fftshift(fft(t_cc_matrix,size(t_cc_matrix,1),1),1);
        figure
        imagesc(Range_bin,f_axis,10*log10(abs(CAF./max(CAF(:)))),[-50 1]); 
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
        

     % Spectrogram 
        int_bins = sum(cc_matrix(2:10,:),1);
        r_bin = 1;
        l_fft = 256;
        pad_factor = 1;
        overlap_factor = 0.99;
        [spect,f] = spectrogram(int_bins,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,passive.seg_s,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1Blackman-Harris:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/passive.Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
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
           passive.doppler_window = 'Blackman-Harris';
           passive.dopp_zero_padding = 1;
           passive.dynamic_range = +50;
           passive.range_bin_size = (1/passive.Fs)/passive.range_zero_padding * C/2;   
           passive.max_range_m = passive.max_range*passive.range_bin_size;
     
                [passive.number_cpi,... 
                 passive.pulses_per_cpi,...
                 passive.range_doppler_slices] = rangeDopplerSlice(cc_matrix,passive.cpi,passive.PRF,...
                                                                   passive.cpi_overlap,...
                                                                   passive.dopp_zero_padding,...
                                                                   passive.doppler_window);  
                
                 [~,~,passive.self_ambg_slices] = rangeDopplerSlice(self_ambg_matrix,passive.cpi,passive.PRF,...
                                                                   passive.cpi_overlap,...
                                                                   passive.dopp_zero_padding,...
                                                                   passive.doppler_window);  
            
            passive.cpi_stride = round(passive.pulses_per_cpi*(1-passive.cpi_overlap)); % number of pulses to stride each for next CPI
            passive.velocity_conv = C*(((1/C)/(passive.Fc/C)));
            passive.range_bins = size(cc_matrix,1);
            passive.doppler_bins = passive.pulses_per_cpi*passive.dopp_zero_padding+1;
            passive.doppler_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.doppler_bins);
            passive.doppler_velocity_axis = passive.doppler_axis*passive.velocity_conv;
            passive.range_axis = linspace(0,passive.max_range_m,passive.range_bins);


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
%                  
%% Direct Signal Interference Cancellation
   % set DSI cancellation parameters
      pfa = 0.999;        % subtraction parameter - P must be a positive 
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
                                               threshold,pfa,...
                                               passive.range_axis,passive.doppler_axis);

% close all
  % create video of CLEANed range-Doppler slices
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

  %% CFAR CAF Slices 

  %% Active CFAR
      % Create look up table to convert range-Doppler cell index to range
      % and Doppler estimate.
        active.range_lookup = ones(size(active.range_doppler_slices{1})).*active.range_axis;
        active.doppler_lookup = ones(size(active.range_doppler_slices{1})).*transpose(active.doppler_velocity_axis);

%         % Create the CFAR detector.
            pfa = 1e-12;
%             active.detector = mod_CFARDetector2D('TrainingBandSize',[4,4], ...
%                                              'GuardBandSize',[3,3],...
%                                              'ThresholdFactor','Auto',...
%                                              'Method','CA', ...
%                                              'ProbabilityFalseAlarm',pfa,...
%                                              'ThresholdOutputPort',true);
        % Create CUT Matrix
            active.full_cutidx = createCutMatrix(active.detector,active.range_bins,size(active.doppler_axis,2));                            
% 
% %         % loop through matrix and operate CFAR to get range-Doppler CFAR results   
%              [active.cfar.range_doppler_slices] = cfarSlices(active.range_doppler_slices,active.detector,active.full_cutidx);
%        
        % Create the CFAR detector.
            active.detector = mod_CFARDetector2D('TrainingBandSize',[4,4], ...
                                             'GuardBandSize',[3,3],...
                                             'ThresholdFactor','Auto',...
                                             'Method','CA', ...
                                             'ProbabilityFalseAlarm',pfa,...
                                             'ThresholdOutputPort',true,...
                                             'NoisePowerOutputPort',true,...
                                             'OutputFormat','Detection index');             
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

%        % create video of CFAR data
%              video_name = exp_dir + "CFAR_active_range-Doppler_Exp_" + Experiment_ID + ".avi";
%              video_title = "CFAR Active Radar Capture";
%              max_range = +150;
%              max_doppler = +inf;
%              frame_rate = 1/(capture_duration/active.number_cpi);   
%        % create video of CLEANed range-Doppler slices
%              createCFARVideo(active.cfar.range_doppler_slices,frame_rate,...
%                              active.range_axis,max_range,...
%                             -active.doppler_axis,max_doppler,...
%                             video_name,video_title);
%     
%% Passive CFAR
%         % Create the CFAR detector.
%             passive.detector = mod_CFARDetector2D('TrainingBandSize',[4,4], ...
%                                              'GuardBandSize',[3,3],...
%                                              'ThresholdFactor','Auto',...
%                                              'Method','CA', ...
%                                              'ProbabilityFalseAlarm',pfa,...
%                                              'ThresholdOutputPort',true);
        % Create CUT Matrix
            passive.full_cutidx = createCutMatrix(passive.detector,passive.range_bins,size(passive.doppler_axis,2));                            

%         % loop through matrix and operate CFAR to get range-Doppler CFAR results   
%              [passive.cfar.range_doppler_slices] = cfarSlices(passive.CLEANed_range_doppler_slices,passive.detector,passive.full_cutidx);
%        
        % Create the CFAR detector.
            passive.detector = mod_CFARDetector2D('TrainingBandSize',[4,4], ...
                                             'GuardBandSize',[3,3],...
                                             'ThresholdFactor','Auto',...
                                             'Method','CA', ...
                                             'ProbabilityFalseAlarm',pfa,...
                                             'ThresholdOutputPort',true,...
                                             'NoisePowerOutputPort',true,...
                                             'OutputFormat','Detection index');             
         % loop through matrix and operate CFAR to get detections, threshold and noise power estimate    
            [passive.cfar.detections,...
             passive.cfar.threshold,...
             passive.cfar.noise_estimate] = cfarSlicesDetections(passive.CLEANed_range_doppler_slices,passive.detector,...
                                                passive.full_cutidx,passive.range_axis,-passive.doppler_velocity_axis);
         
         % plot detections in range and Doppler
            figure
            plot(interp_flight_path_cpi,2*-interp_sog_cpi,'black-',LineWidth=2)
            hold on
            for i=1:active.number_cpi
                plot(active.cfar.detections{i}(1,:),active.cfar.detections{i}(2,:),'r*')
                hold on
                plot(passive.cfar.detections{i}(1,:)+10,passive.cfar.detections{i}(2,:),'b*')
                hold on
            end
            grid on
            ylabel('Doppler Velocity (m/s)')
            xlabel('Range (m)')  
            xlim([-inf 100])
            ylim([-20 20])
            legend('Active Radar Detections','Passive Radar Detections')


%        % create video of CFAR data
%              video_name = exp_dir + "CFAR_passvie_range-Doppler_Exp_" + Experiment_ID + ".avi";
%              video_title = "CFAR Active Radar Capture";
%              max_range = 200;
%              max_doppler = 200;
%              frame_rate = 1/(capture_duration/passive.number_cpi);    
%  
% % create video of CFARed CLEANed range-Doppler slices
%              createCFARVideo(passive.cfar.range_doppler_slices,frame_rate,...
%                              passive.range_axis,max_range,...
%                             -passive.doppler_axis,max_doppler,...
%                             video_name,video_title);


    

%% Hybrid Video 
 
   video_name = exp_dir + "hybrid_range-Doppler_Exp_" + Experiment_ID + ".avi";     
   dynamic_range = 50;
   range_limit = 200;
   doppler_limit = 20;
            createVideos(active.range_doppler_slices,...
                         passive.CLEANed_range_doppler_slices,...
                         frame_rate,...
                         active.range_axis,...
                         passive.range_axis,...
                         -active.doppler_velocity_axis,...
                         passive.doppler_velocity_axis,...
                         range_limit,...
                         doppler_limit,...
                         dynamic_range,...
                         video_name,video_title);
                     
%% Hybrid Range and Doppler Plots 

number_cpi = size(passive.inter_range_doppler_slices,2);
passive_range_detections = zeros(number_cpi,1);
active_range_detections = zeros(number_cpi,1);
passive_doppler_detections = zeros(number_cpi,1);
active_doppler_detections = zeros(number_cpi,1);
                     
passive_range_lookup = ones(size(passive.inter_range_doppler_slices {1})).*passive.interp_range_axis;
passive_doppler_lookup = ones(size(passive.inter_range_doppler_slices{1})).*transpose(passive.interp_doppler_velocity_axis);

active_range_lookup = ones(size(active.range_doppler_slices{1})).*active.range_axis;

active_doppler_lookup = ones(size(active.range_doppler_slices{1})).*transpose(active.doppler_velocity_axis);

for i=1:number_cpi   
   i
   [mag, index] = max(passive.inter_range_doppler_slices{i}(:));
   passive_range_detections(i) = passive_range_lookup(index);
   passive_doppler_detections(i) = passive_doppler_lookup(index);
   
   
   [mag, index] = max(active.range_doppler_slices{i}(:));
   active_range_detections(i) = active_range_lookup(index);
   active_doppler_detections(i) = active_doppler_lookup(index);                  
end

time_axis = linspace(0,capture_duration,number_cpi);

fig = figure
pfa = plot(time_axis,passive_range_detections+10);
hold on
a = plot(time_axis,active_range_detections-2);
pfa.Marker = '*';
a.Marker = '*';
% plot(0:capture_duration,gps_range,'r-')
grid on
grid minor
ylabel('Range (m)')
xlabel('Time (s)')
% ylim([-inf 50])
legend('Passive Radar Peak Return','Active Radar Peak Return','GPS Ground Truth')
fig_name = exp_dir + "Hybrid Data RT - " + Experiment_ID + ".jpg";
saveas(fig,fig_name,'jpeg')
saveas(fig,fig_name)


fig = figure
plot(time_axis,-passive_doppler_detections,'-bl*');
hold on
a = plot(time_axis,active_doppler_detections,'-r*');
plot(0:capture_duration,2*-gps_sog,'black-')
grid on
grid minor
% ylim([-8 8])
ylabel('Velocity (m/s)')
xlabel('Time (s)')
legend('Passive Radar Peak Return','Active Radar Peak Return','GPS Ground Truth')
fig_name = exp_dir + "Hybrid Data DT - " + Experiment_ID + ".jpg";
saveas(fig,fig_name,'jpeg')
saveas(fig,fig_name)


% end


%% SNR of peak return 

passive_range_detections = zeros(active.number_cpi,1);
passive_detection_snr = zeros(active.number_cpi,1);
active_range_detections = zeros(active.number_cpi,1);
active_detection_snr = zeros(active.number_cpi,1);
for i=1:active.number_cpi
   i
   [target_mag, index] = max(passive.CLEANed_range_doppler_slices{i}(:));
   passive_range_detections(i) = passive_range_lookup(index)+10
   passive_target_power = 10*log10(abs(target_mag))
   passive_noise_estimate = 10*log10(mean(abs(passive.CLEANed_range_doppler_slices{i}(:,45:50)),'all'))
   passive_detection_snr(i) = passive_target_power - passive_noise_estimate

   [target_mag, index] = max(active.range_doppler_slices{i}(:));
   active_range_detections(i) = active_range_lookup(index)
   active_target_power = 10*log10(abs(target_mag))
   active_noise_estimate = 10*log10(mean(abs(active.range_doppler_slices{i}(:,150:170)),'all'))
   active_detection_snr(i) = active_target_power - active_noise_estimate

end  
t = linspace(0,15,size(passive_detection_snr,1));
fig = figure    
plot(passive_range_detections,passive_detection_snr,'*r')    
hold on
plot(active_range_detections,active_detection_snr,'*b')    
legend('SNR Passive Radar Peak Return','SNR Active Radar Peak Return')
grid on
grid minor
ylabel('SNR (dB)')
xlabel('Range (m)')
% xlim([-inf 100])
fig_name = exp_dir + "Passive SNR of Peak Detection - CPI = " + passive.cpi + " s.jpg";
saveas(fig,fig_name,'jpeg')
saveas(fig,fig_name)

%% Plot GPS Ground Truth Data 
    gps_file_time_offset = 280-14; % offset between start time of radar and gps recording
    gps_file_range_offset = 8;%8 
%     gps_file_time_offset = 3120-11; % offset between start time of radar and gps recording
%     gps_file_range_offset = 8;
    gps_range = geode1.target_range(gps_file_time_offset:gps_file_time_offset+capture_duration) + gps_file_range_offset;
    gps_sog = geode1.SoG(gps_file_time_offset:gps_file_time_offset+capture_duration);
%     figure
%     subplot(2,1,1)
%     plot(0:15,gps_range)
%     ylabel('Range (m)')
%     xlabel('Time (s)')
%     grid on 
%     subplot(2,1,2)
%     plot(0:15,gps_sog)
%     ylabel('Speed over Ground (m/s)')
%     xlabel('Time (s)')
%     grid on 
% Interpolate to number of pulses
    old_sample_points = 0:capture_duration;
    active_pulses_sample_points = 1:active.number_pulses;
    interp_flight_path = interp1(old_sample_points,gps_range,active_pulses_sample_points/1000,'linear');
    %     figure
    %     plot(interp_drone_flight_path)
    %     xlabel('Range (m)')
    %     ylabel('Time (s)')
% Interpolate to number of CPI
       cpi_sample_points = linspace(0,capture_duration,active.number_cpi);
       interp_flight_path_cpi = interp1(old_sample_points,gps_range,cpi_sample_points,'linear');
       interp_sog_cpi = interp1(old_sample_points,gps_sog,cpi_sample_points,'linear');

               
% Range Limited Plot 
    fig = figure
    imagesc(active.range_axis,active.time_axis,MTI_RTI_plot,[-30,0]);               
    hold on 
    plot(interp_flight_path,active.time_axis,'-r', 'LineWidth', 2);
    grid on            
    c = colorbar
    c.Label.String='Norm Power (dB)'
    ylabel('Time (Sec)')
    xlabel('Range (m)')            
    xlim([0 200])
    legend('Target GPS Ground Truth')

    fig_name = exp_dir + "Active_radar_capture_vs_ground_truth_" + Experiment_ID;
    saveas(fig,fig_name,'jpeg')
    saveas(fig,fig_name)
        
% Range Limited Plot 
    fig = figure
    imagesc(active.range_axis,active.time_axis,MTI_RTI_plot,[-50,0]);               
    hold on 
    plot(interp_drone_flight_path,time_axis,'rd', 'MarkerSize', 2);
    grid on            
    colorbar
    ylabel('Time (Sec)')
    xlabel('Range Bin')            
    xlim([0 100])
    fig_name = exp_dir + "MTI_active_radar_capture_vs_ground_truth_" + Experiment_ID;
    saveas(fig,fig_name,'jpeg')
    saveas(fig,fig_name)

% plot detections and ground truth in range and Doppler
        figure
        plot(interp_flight_path_cpi,2*-interp_sog_cpi,'black-','LineWidth',2)
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

%% Pd
% compute the active sensor probability of detection
active.cfar.range_error_lim = 2.5; 
active.cfar.sog_error_lim = 2;
active.rangeBias = -2;
[active.cfar.detected, active.cfar.pd] = computePD(active.cfar.detections,interp_flight_path_cpi,2*-interp_sog_cpi,active.cfar.range_error_lim,active.cfar.sog_error_lim,active.rangeBias);
figure; plot(active.cfar.detected); ylim([0 2]);title('Active Detections');active.cfar.pd

% compute the passive sensor probability of detection
passive.cfar.range_error_lim = 7.5; 
passive.cfar.sog_error_lim = 1;
passive.rangeBias = 8;
[passive.cfar.detected, passive.cfar.pd] = computePD(passive.cfar.detections,interp_flight_path_cpi,2*-interp_sog_cpi,passive.cfar.range_error_lim,passive.cfar.sog_error_lim,passive.rangeBias);
figure; plot(passive.cfar.detected); ylim([0 2]); title('Passive Detections'); passive.cfar.pd
 
% OR fused detection result
fused = active.cfar.detected + passive.cfar.detected;
fused = ceil(fused/2);
or_fused_pd = (sum(fused)/active.number_cpi)*100

% AND fused detection result
fused = active.cfar.detected + passive.cfar.detected;
fused = floor(fused - 0.5);
and_fused_pd = (sum(fused)/active.number_cpi)*100


figure; plot(fused); ylim([0 2]);title('Fused Detections');






% passive.xcorr_len = 0.5; 
%     passive.nSurfs = 60;
%     passive.doppler_range = 100;
%     passive.doppler_resolution = 1;
%     passive.range_bins = 50;
% 
%     passive.range_doppler_slices = passive_process(ref_channel,sur_channel,passive.xcorr_len,passive.nSurfs,passive.doppler_range,passive.doppler_resolution,passive.Fs,passive.range_bins,exp_dir,passive.range_zero_padding,passive.td_corr)
%          
%              


g = 30;
[nAdbins nArdbins] = size(active.range_doppler_slices{g});
[nPdbins nPrdbins] = size(passive.inter_range_doppler_slices{g});
active_data = flipud(active.range_doppler_slices{g}(nAdbins/2-50:nAdbins/2+50,1:100));
passive_data = passive.inter_range_doppler_slices{g}(nPdbins/2-50:nPdbins/2+50,1:100);
active_data = abs(active_data./max(active_data(:)));
passive_data = abs(passive_data./max(passive_data(:)));
fused_data = active_data + passive_data;

data = active_data;
frame =  10*log(abs(data./max(data(:))));
figure
imagesc(frame, [-50 0])

data = passive_data;
frame =  10*log(abs(data./max(data(:))));
figure
imagesc(frame, [-50 0])

data = fused_data;
frame =  10*log(abs(data./max(data(:))));
figure
imagesc(frame, [-50 0])



fused = active.range_doppler_slices{1}(1:100,125:226) ;+ size(passive.inter_range_doppler_slices{1});