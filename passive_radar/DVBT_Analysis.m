%% This file is for analysis of DVB-T radar captures 





%% How to use analysis scripts 
% 1) load the 'experimental configuration.mat' file
% 2) enter correct location of experiments as save_directory1

% location of the raw data
save_directory1 = "/media/piers/T7/06_07_2022_farm/hybrid_radar/"
exp_dir = save_directory1 + Experiment_ID + '/';

% clear all
close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/',...
        '~/repos/bladeRAD/generic_scripts/matlab/CFAR') % path to generic functions

%% Capture parameters 
passive.max_range = 100 % max number of range bins to cross correlate
passive.bin_size = C*1/(passive.Bw*2)

 %% Passive Processing
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
         passive.seg_s = 1000; % number of segments per second - analagos to PRF.
         passive.seg_percent = 100;  % percentage of segment used for cross coreclation of 
                            % survallance and reference. Will affect SNR dramatically.
         [ref_matrix ,self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.max_range,exp_dir);
         save(exp_dir + 'passive_matrix','cc_matrix')
    % RTI Plot
        RTI_plot= transpose(10*log10(abs(cc_matrix./max(cc_matrix(:)))));
        passive.range_bins = size(cc_matrix,1);  
        passive.range_axis = linspace(0,passive.max_range*passive.bin_size,passive.range_bins);
        time_axis = linspace(0,capture_duration,size(cc_matrix,2));
        figure
        fig = imagesc(passive.range_axis,time_axis,RTI_plot,[-inf,0]);
            % xlim([1 20])
            %ylim([0 0.0005])
            grid on            
            c = colorbar;
            c.Label.String='Norm Power (dB)';
            ylabel('Time (Sec)')
            xlabel('Range (m)')   
            fig_title = "Passive RTI - " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "Passive RTI_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)

      % CAF of entire capture
        f_axis = linspace(-passive.seg_s/2,passive.seg_s/2,size(cc_matrix,2));
        t_cc_matrix = transpose(cc_matrix);
        CAF = fftshift(fft(t_cc_matrix,size(t_cc_matrix,1),1),1);
        fig = figure
        imagesc(passive.range_axis,f_axis,10*log10(abs(CAF./max(CAF(:)))),[-50 1]); 
            ylim([-50 50])     
            xlim([-inf 400])
            c = colorbar;
            c.Label.String='Norm Power (dB)';
            ylabel('Doppler Shift (Hz)')
            xlabel('Range (m)')  
            title("CAF for entire capture" + Experiment_ID)
            fig_name = exp_dir + "CAF for entire capture_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)
        

     % Spectrogram 
        r_start = 1;
        r_stop = 100;
        l_fft = 1024;
        pad_factor = 4;
        overlap_factor = 0.99;
        integrated_data = sum(cc_matrix(r_start:r_stop,:));
        [spect,f] = spectrogram(integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,passive.seg_s,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/passive.Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-50 0]);   
            ylim([-50 50])
            colorbar
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "Passive Spectrogram - " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "Passive Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)


%% Proccess Passive data into Range-Doppler Slices
           passive.PRF = passive.seg_s; %seg_s
           passive.cpi = 0.5; % cohernet proccessing interval (s)
           passive.cpi_overlap = 0.8; % overlap between CPIs (watch this - too large will cause slow execution)
           passive.doppler_window = 'hann';
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
            passive.pulses_per_cpi

            passive.cpi_stride = round(passive.pulses_per_cpi*(1-passive.cpi_overlap)); % number of pulses to stride each for next CPI
            passive.velocity_conv = C*(((1/C)/(passive.Fc/C)));
            passive.doppler_bins = passive.pulses_per_cpi*passive.zero_padding+1;
            passive.doppler_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.doppler_bins);
            passive.doppler_velocity_axis = passive.doppler_axis*passive.velocity_conv;

%             % create video of CLEANed range-Doppler slices
%              video_name = exp_dir + "Range_Doppler_Slices" + Experiment_ID + ".avi";
%              %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";       
%              video_title = "Passive Pre-DSI";
%              dynamic_range = passive.dynamic_range;
%              max_range = 500;
%              max_doppler = 100;
%              frame_rate = 1/(capture_duration/passive.number_cpi);    
%              createVideo(passive.range_doppler_slices,frame_rate,...
%                          passive.range_axis,max_range,...
%                          passive.doppler_axis,max_doppler,...
%                          dynamic_range,video_name,video_title);
%         
                 
%% Direct Signal Interference Cancellation
   % set DSI cancellation parameters
      p = 0.99;        % subtraction parameter - P must be a positive 
                        % integer less than one to avoid    unwanted discontinuities
                        % arising from zero values in the range–Doppler surface.
     threshold = 0.005; % cutoff threshold parameter
     max_iterations = 500; % maximum number of itterations DSI is CLEANed from CAF slice
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


%   % create video of CLEANed range-Doppler slices
%      video_name = exp_dir + "CLEANed_range-Doppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";
%      %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";       
%      video_title = "CLEANed Passive Radar Capture";
%      dynamic_range = 30;
%      max_range = 500;
%      max_doppler = 100;
%      frame_rate = 1/(capture_duration/passive.number_cpi);    
%      createVideo(passive.CLEANed_range_doppler_slices,frame_rate,...
%                  passive.range_axis,max_range,...
%                  passive.doppler_axis,max_doppler,...
%                  dynamic_range,video_name,video_title);


%% maximum return detector
number_cpi = size(passive.CLEANed_range_doppler_slices,2);
passive_range_lookup = ones(size(passive.CLEANed_range_doppler_slices{1})).*passive.range_axis;
passive_doppler_lookup = ones(size(passive.CLEANed_range_doppler_slices{1})).*transpose(passive.doppler_velocity_axis);
time_axis = linspace(0,15,number_cpi);


passive_range_max_det = zeros(number_cpi,1);
passive_doppler_max_det = zeros(number_cpi,1);
for i=1:number_cpi   
   i
   [mag, index] = max(passive.CLEANed_range_doppler_slices{i}(:));
   passive_range_max_det(i) = passive_range_lookup(index);
   passive_doppler_max_det(i) = passive_doppler_lookup(index);            
end


fig = figure
p = plot(time_axis,passive_range_max_det);
p.Marker = '*';
% a.Marker = '*';
grid on
grid minor
ylabel('Range (m)')
xlabel('Time (s)')
ylim([-inf 100])
% ylim([-inf 50])
legend('Passive Radar Peak Return')
fig_name = exp_dir + "Passive Peak Detection- CPI = " + passive.cpi + " s.jpg";
saveas(fig,fig_name,'jpeg')
saveas(fig,fig_name)

%% SNR of peak return 

passive_range_detections = zeros(number_cpi-1,1);
passive_detection_snr = zeros(number_cpi-1,1);
for i=1:number_cpi-1
   i
   [target_mag, index] = max(passive.CLEANed_range_doppler_slices{i}(:));
   passive_range_detections(i) = passive_range_lookup(index)
   passive_target_power = 10*log10(abs(target_mag))
   passive_noise_estimate = 10*log10(mean(abs(passive.CLEANed_range_doppler_slices{i}(:,70:101)),'all'))
   passive_detection_snr(i) = passive_target_power - passive_noise_estimate
end  
t = linspace(0,15,size(passive_detection_snr,1));
fig = figure    
plot(passive_range_detections,passive_detection_snr,'*r')    
grid on
grid minor
ylabel('SNR (dB)')
xlabel('Range (m)')
xlim([-inf 100])
fig_name = exp_dir + "Passive SNR of Peak Detection - CPI = " + passive.cpi + " s.jpg";
saveas(fig,fig_name,'jpeg')
saveas(fig,fig_name)


%% CFAR CAF Slices 
        % Create the CFAR detector.
            p = 1e-6;
            passive_detector = mod_CFARDetector2D('TrainingBandSize',[8,8], ...
                                             'GuardBandSize',[6,6],...
                                             'ThresholdFactor','Auto',...
                                             'Method','CA', ...
                                             'ProbabilityFalseAlarm',p,...
                                             'ThresholdOutputPort',true);
        % Create CUT Matrix
            passive.full_cutidx = createCutMatrix(passive_detector,passive.range_bins,passive.doppler_bins);                             

        %loop through matrix and operate CFAR     
            [passive.cfar_range_doppler_slices] = cfarSlices(passive.CLEANed_range_doppler_slices,passive_detector,passive.full_cutidx);
   
       % create video of CFAR data
             video_name = exp_dir + "CFAR_passive_range-Doppler_Exp_" + Experiment_ID + ".avi";
             video_title = "CFAR Passive Radar Capture";
             max_range = 200;
             max_doppler = 200;
             frame_rate = 1/(capture_duration/passive.number_cpi);    

             createCFARVideo(passive.cfar_range_doppler_slices,frame_rate,...
                             passive.range_axis,max_range,...
                            -passive.doppler_axis,max_doppler,...
                            video_name,video_title);

