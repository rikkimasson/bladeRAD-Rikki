%% Fuse active and passive range Doppler images 

%% Interpolate passive radar data
    range_factor = passive.range_bin_size/active.range_bin_size; % Interolation factor in range
    doppler_factor = passive.velocity_conv/active.velocity_conv; % Interpolation factor in doppler 
    
    passive.interp_range_doppler_slices = interpolateRangeDoppler(passive.CLEANed_range_doppler_slices,range_factor,doppler_factor);
  
    [passive.interp_doppler_bins passive.interp_range_bins] = size(passive.interp_range_doppler_slices{1});
    passive.interp_doppler_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.interp_doppler_bins);
    passive.interp_doppler_velocity_axis = passive.interp_doppler_axis*passive.velocity_conv;
    passive.interp_range_axis = linspace(0,passive.max_range_m,passive.interp_range_bins);


    % create video of  range-Doppler slices
     video_name = exp_dir + "Interp_Range_Doppler_Slices" + Experiment_ID + ".avi";
     video_title = "";
     dynamic_range = +80;
     max_range = 200;
     max_doppler = 100;
     frame_rate = 1/(capture_duration/passive.number_cpi);    
     createVideo(passive.interp_range_doppler_slices,frame_rate,...
                 passive.range_axis,max_range,...
                 passive.doppler_axis,max_doppler,...
                 dynamic_range,video_name,video_title);


%% Scale range Doppler images and fuse
passive_range_slide = 0;
active_range_slide  = 2; %2
[nAdbins nArbins] = size(active.range_doppler_slices{1});
[nPdbins nPrbins] = size(passive.interp_range_doppler_slices{1});
hybrid.no_range_bins = min([nArbins nPrbins])-max([active_range_slide passive_range_slide]);
hybrid.no_doppler_bins = min([nAdbins nPdbins]);
hybrid.no_doppler_bins_2 = hybrid.no_doppler_bins/2-1;
hybrid.range_doppler_slices = createArrays(size(active.range_doppler_slices,2), [hybrid.no_range_bins hybrid.no_doppler_bins]);
hybrid.range_axis = active.range_axis(1:hybrid.no_range_bins);
hybrid.doppler_velocity_axis = active.doppler_velocity_axis(nAdbins/2-hybrid.no_doppler_bins_2:nAdbins/2+hybrid.no_doppler_bins_2+1);


wv = 'db2';
lv = 5;
active_data = createArrays(size(active.range_doppler_slices,2), [hybrid.no_range_bins hybrid.no_doppler_bins]);
passive_data = createArrays(size(active.range_doppler_slices,2), [hybrid.no_range_bins hybrid.no_doppler_bins]);

for i=1:size(active.range_doppler_slices,2)
    active_data{i} = flipud(active.range_doppler_slices{i}(nAdbins/2-hybrid.no_doppler_bins_2:nAdbins/2+hybrid.no_doppler_bins_2+1,1+active_range_slide:hybrid.no_range_bins+active_range_slide));
    passive_data{i} = passive.interp_range_doppler_slices{i}(nPdbins/2-hybrid.no_doppler_bins_2:nPdbins/2+hybrid.no_doppler_bins_2+1,1+passive_range_slide:hybrid.no_range_bins+passive_range_slide);
    active_data1 = (active_data{i}./max(active_data{i}(:)));
    passive_data1 = (passive_data{i}./max(passive_data{i}(:)));
%     hybrid.range_doppler_slices{i} = wfusmat(abs(active_data{i}),abs(passive_data{i}),'mean');
    [aA,aH,aV,aD] = dwt2(abs(active_data{i}),'sym5','mode','per');
    [pA,pH,pV,pD] = dwt2(abs(passive_data{i}),'sym5','mode','per');
    
%     hybrid.range_doppler_slices{i} = idwt2((aA),(aH+pH),(aV+pV),(aD+pD),'db4')

    hybrid.range_doppler_slices{i} = abs(active_data1) + abs(passive_data1);

end 

% % create video of  range-Doppler slices
% %      video_name = exp_dir + "hybrid.no_range_bins_Doppler_Slices" + Experiment_ID + ".avi";
% %      video_title = "";
% %      dynamic_range = +90;
% %      max_range = +inf;
% %      max_doppler = 60;
% %      frame_rate = 1/(capture_duration/passive.number_cpi);    
% %      createVideo(hybrid.range_doppler_slices,frame_rate,...
% %                  hybrid.range_axis,max_range,...
% %                  hybrid.doppler_velocity_axis,max_doppler,...
% %                  dynamic_range,video_name,video_title);

  % create video of  range-Doppler slices
    video_name = exp_dir + "Three images" + Experiment_ID + ".avi";
    video_title = "";
    dynamic_range = +100;
    max_range = +inf;
    max_doppler = 60;
    createVideo3(active_data, ...
                passive_data,hybrid.range_doppler_slices, ...
                frame_rate,hybrid.range_axis,hybrid.doppler_velocity_axis, ...
                max_range,max_doppler,dynamic_range,video_name,video_title)



%% SNR of peak return 
% passive_range_lookup = ones(size(passive.CLEANed_range_doppler_slices{2})).*passive.range_axis;
% passive_doppler_lookup = ones(size(passive.CLEANed_range_doppler_slices{2})).*transpose(passive.doppler_velocity_axis);

passive_range_lookup = ones(size(passive.interp_range_doppler_slices {2})).*passive.interp_range_axis;
passive_doppler_lookup = ones(size(passive.interp_range_doppler_slices{2})).*transpose(passive.interp_doppler_velocity_axis);

active_range_lookup = ones(size(active.range_doppler_slices{2})).*active.range_axis;
active_doppler_lookup = ones(size(active.range_doppler_slices{2})).*transpose(active.doppler_velocity_axis);

hybrid_range_lookup = ones(size(hybrid.range_doppler_slices{2})).*hybrid.range_axis;
hybrid_doppler_lookup = ones(size(hybrid.range_doppler_slices{2})).*transpose(hybrid.doppler_velocity_axis);


passive_range_detections = zeros(active.number_cpi,1);
passive_detection_snr = zeros(active.number_cpi,1);
active_range_detections = zeros(active.number_cpi,1);
active_detection_snr = zeros(active.number_cpi,1);
hybrid_range_detections = zeros(active.number_cpi,1);
hybrid_detection_snr = zeros(active.number_cpi,1);
for i=2:active.number_cpi
   i
%    [target_mag, index] = max(passive.CLEANed_range_doppler_slices{i}(:));
%    passive_range_detections(i) = passive_range_lookup(index)+10
%    passive_target_power = 10*log10(abs(target_mag))
%    passive_noise_estimate = 10*log10(mean(abs(passive.CLEANed_range_doppler_slices{i}(:,40:51)),'all'))
%    passive_detection_snr(i) = passive_target_power - passive_noise_estimate

   [target_mag, index] = max(passive.interp_range_doppler_slices{i}(:));
   passive_range_detections(i) = passive_range_lookup(index);
   passive_target_power = 10*log10(abs(target_mag));
   passive_noise_estimate = 10*log10(mean(abs(passive.interp_range_doppler_slices{i}(350:500,110:120)),'all'));   
   passive_detection_snr(i) = passive_target_power - passive_noise_estimate;


   [target_mag, index] = max(active.range_doppler_slices{i}(:));
   active_range_detections(i) = active_range_lookup(index);
   active_target_power = 10*log10(abs(target_mag));
   active_noise_estimate = 10*log10(mean(abs(active.range_doppler_slices{i}(350:500,80:100)),'all'));
   active_detection_snr(i) = active_target_power - active_noise_estimate;

   [target_mag, index] = max(hybrid.range_doppler_slices{i}(:));
   hybrid_range_detections(i) = hybrid_range_lookup(index);
   hybrid_target_power = 10*log10(abs(target_mag));
   hybrid_noise_estimate = 10*log10(mean(abs(hybrid.range_doppler_slices{i}(:,78:98)),'all'));
   hybrid_detection_snr(i) = hybrid_target_power - hybrid_noise_estimate;

% figure
% imagesc(10*log10(abs(active.range_doppler_slices{58}))); colorbar;colormap jet;


end  
t = linspace(0,15,size(passive_detection_snr,1));
fig = figure    
plot(passive_range_detections,passive_detection_snr,'*r')    
hold on
plot(active_range_detections-10,active_detection_snr,'*b')    
hold on
plot(hybrid_range_detections,hybrid_detection_snr,'*magenta')   
legend('SNR Passive Radar Peak Return','SNR Active Radar Peak Return','SNR of Hybrid Peak Return')
grid on
grid minor
ylabel('SNR (dB)')
xlabel('Range (m)')
% xlim([-inf 100])
fig_name = exp_dir + "Passive SNR of Peak Detection - CPI = " + passive.cpi + " s.jpg";
saveas(fig,fig_name,'jpeg')
saveas(fig,fig_name)




  %% CFAR Interpolated Passive Data
    pfa = 1e-12;     
    % There are two CFAR detectors        
               % 1. Provides the detections cell array 'active.cfar.detections', 
               %    each cell is vector per CPI of detections in range and doppler.
               % 2. Provides cell array of binary valued range-Doppler surfaces 
               %    with detections as 1.
    
         % 1. Array of vectors with detection in range and Doppler   
            % Create the CFAR detector.
                hybrid.detector = mod_CFARDetector2D('TrainingBandSize',[8,8], ...
                                                 'GuardBandSize',[6,6],...
                                                 'ThresholdFactor','Auto',...
                                                 'Method','CA', ...
                                                 'ProbabilityFalseAlarm',pfa,...
                                                 'ThresholdOutputPort',true,...
                                                 'NoisePowerOutputPort',true,...
                                                 'OutputFormat','Detection index');             
           % Create CUT Matrix
                hybrid.full_cutidx = createCutMatrix(hybrid.detector,hybrid.no_range_bins,hybrid.no_doppler_bins);                            
     
           % Loop through matrix and operate CFAR to get detections, threshold and noise power estimate    
                [hybrid.cfar.detections,...
                 hybrid.cfar.threshold,...
                 hybrid.cfar.noise_estimate] = cfarSlicesDetections(hybrid.range_doppler_slices,hybrid.detector,...
                                                     hybrid.full_cutidx,hybrid.range_axis,-hybrid.doppler_velocity_axis);
             
                     
           % Plot Hybrid detections in range and Doppler
                figure
                plot(gt.interp_range-10,2*-gt.interp_sog,'black-',LineWidth=2);hold on;      
                for i=1:active.number_cpi
                    plot(hybrid.cfar.detections{i}(1,:),hybrid.cfar.detections{i}(2,:),'m*');hold on;        
                end
                grid on; ylabel('Doppler Velocity (m/s)'); xlabel('Range (m)'); legend('Ground Truth','Hybrid Radar Detections');
                xlim([0 100]); ylim([-20 20]);
                
          % 2. Array of vectors with detection in range and Doppler
            % Create the CFAR detector.
                hybrid.detector = mod_CFARDetector2D('TrainingBandSize',[8,8], ...
                                                 'GuardBandSize',[6,6],...
                                                 'ThresholdFactor','Auto',...
                                                 'Method','CA', ...
                                                 'ProbabilityFalseAlarm',pfa,...
                                                 'ThresholdOutputPort',true);
            
            % loop through matrix and operate CFAR to get range-Doppler CFAR results   
                 [hybrid.cfar.range_doppler_slices] = cfarSlices(hybrid.range_doppler_slices,hybrid.detector,hybrid.full_cutidx);
           
       
        % create video of CFAR data
             video_name = exp_dir + "CFAR_passvie_range-Doppler_Exp_" + Experiment_ID + ".avi";
             video_title = "CFAR Active Radar Capture";
             max_range = 200;
             max_doppler = 200;
             frame_rate = 1/(capture_duration/passive.number_cpi);    
             figure;
             createCFARVideo(hybrid.cfar.range_doppler_slices,frame_rate,...
                             hybrid.range_axis,max_range,...
                             hybrid.doppler_axis,max_doppler,...
                             video_name,video_title);

%% Calculate Scenario Pd
    % Currently 'computePD' function is setup just to check if there is a 
    % detection at the correct range, not range and Doppler.
    
    % compute the active sensor probability of detection
        active.cfar.range_error_lim = 2.5; 
        active.cfar.sog_error_lim = 1;
        active.rangeBias = -2;
        [active.cfar.detected, active.cfar.pd] = computePD(active.cfar.detections,gt.interp_range,2*-gt.interp_sog,active.cfar.range_error_lim,active.cfar.sog_error_lim,active.rangeBias);
%         figure; plot(active.cfar.detected); ylim([0 2]);title('Active Detections');active.cfar.pd
    
    % compute the passive sensor probability of detection
        passive.cfar.range_error_lim = 7.5; 
        passive.cfar.sog_error_lim = 1;
        passive.rangeBias = hybrid.range_shift ;
        [passive.cfar.detected, passive.cfar.pd] = computePD(passive.cfar.detections,gt.interp_range,2*-gt.interp_sog,passive.cfar.range_error_lim,passive.cfar.sog_error_lim,passive.rangeBias);
%         figure; plot(passive.cfar.detected); ylim([0 2]); title('Passive Detections'); passive.cfar.pd

    % decentalised detection result
        fused = active.cfar.detected + passive.cfar.detected;
        fused = ceil(fused/2);
        fused_pd = (sum(fused)/active.number_cpi)*100
%         figure; plot(fused); ylim([0 2]);title('De-centralised Detections');
    
    % compute the Centralised sensor probability of detection
        hybrid.cfar.range_error_lim = 5; 
        hybrid.cfar.sog_error_lim = 2;
        hybrid.rangeBias = hybrid.range_shift;
        [hybrid.cfar.detected, hybrid.cfar.pd] = computePD(hybrid.cfar.detections,gt.interp_range,2*-gt.interp_sog,hybrid.cfar.range_error_lim,hybrid.cfar.sog_error_lim,hybrid.rangeBias);
%         figure; plot(hybrid.cfar.detected); ylim([0 2]); title('Centralised Detections'); hybrid.cfar.pd

        figure
        subplot(4,1,1)
        plot(active.cfar.detected,'r-*'); grid minor; 
        title('Active Detections');
        subplot(4,1,2)
        plot(passive.cfar.detected,'b-*'); grid minor; 
        title('Passive Detections');
        subplot(4,1,3)
        plot(fused,'m-*'); hold on;ylim([0 1])
        title('De-centralised Detections'); grid minor; 
        subplot(4,1,4)
        plot(hybrid.cfar.detected,'m-*'); hold on; grid minor; 
        title('Centralised Detections');

    

%% Plot range and Doppler of peak cell in active and passive range-Doppler surfaces.
  % Create empty vectors for detections 
    passive_range_detections = zeros(active.number_cpi,1);
    active_range_detections = zeros(active.number_cpi,1);
    passive_doppler_detections = zeros(active.number_cpi,1);
    active_doppler_detections = zeros(active.number_cpi,1);
    hybrid_doppler_detections = zeros(active.number_cpi,1);
    hybrid_doppler_detections = zeros(active.number_cpi,1);
  % Create lookup matrix to find range and doppler of detection on matrix               
    passive_range_lookup = ones(size(passive.interp_range_doppler_slices{1})).*passive.interp_range_axis;
    passive_doppler_lookup = ones(size(passive.interp_range_doppler_slices{1})).*transpose(passive.interp_doppler_velocity_axis);
    active_range_lookup = ones(size(active.range_doppler_slices{1})).*active.range_axis;
    active_doppler_lookup = ones(size(active.range_doppler_slices{1})).*transpose(active.doppler_velocity_axis);
    hybrid_range_lookup = ones(size(hybrid.range_doppler_slices{1})).*hybrid.range_axis;
    hybrid_doppler_lookup = ones(size(hybrid.range_doppler_slices{1})).*transpose(hybrid.doppler_velocity_axis);
  
   % Find peak cell on range doppler surface for each CPI
    for j=1:active.number_cpi   
       [mag, index] = max(passive.interp_range_doppler_slices{j}(:));
       passive_range_detections(j) = passive_range_lookup(index);
       passive_doppler_detections(j) = passive_doppler_lookup(index);
    
       [mag, index] = max(active.range_doppler_slices{j}(:));
       active_range_detections(j) = active_range_lookup(index);
       active_doppler_detections(j) = active_doppler_lookup(index);         

       [mag, index] = max(hybrid.range_doppler_slices{j}(:));
       hybrid_range_detections(j) = hybrid_range_lookup(index);
       hybrid_doppler_detections(j) = hybrid_doppler_lookup(index);  
    end

    passive.range_shift = 2 * passive.interp_range_axis(2);
    hybrid.range_shift = 2 * passive.interp_range_axis(2);
    % Plot peak detection in range over time   
        time_axis = linspace(0,capture_duration,active.number_cpi);
        fig = figure
        pa = plot(time_axis,passive_range_detections+passive.range_shift); hold on;
        a = plot(time_axis,active_range_detections);hold on;
        c = plot(time_axis,hybrid_range_detections+hybrid.range_shift,'m-*');hold on;
        pfa.Marker = '*'; a.Marker = '*';
        plot(0:capture_duration,gt.range,'black-')
        grid on; grid minor; ylabel('Range (m)'); xlabel('Time (s)');
        % ylim([-inf 50])
        legend('Passive Radar Peak Return','Active Radar Peak Return','Centralised Peak Return','GPS Ground Truth')
        fig_name = exp_dir + "Hybrid Data RT - " + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg'); saveas(fig,fig_name);
    
        std((active_range_detections(7:end)-(passive_range_detections(7:end)+passive.range_shift)))

    % Plot peak detection in velocity over time   
        fig = figure
        plot(time_axis,-passive_doppler_detections,'-bl*'); hold on;
        a = plot(time_axis,active_doppler_detections,'-r*'); hold on;
        c = plot(time_axis,-hybrid_doppler_detections,'m-*');hold on;
        plot(0:capture_duration,2*-gt.sog,'black-')
        grid on; grid minor;
        % ylim([-8 8])
        ylabel('Velocity (m/s)'); xlabel('Time (s)');
        legend('Passive Radar Peak Return','Active Radar Peak Return','Centralised Peak Return','GPS Ground Truth')
        fig_name = exp_dir + "Hybrid Data DT - " + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg'); saveas(fig,fig_name);
    
        mean((active_doppler_detections(7:end)+(passive_doppler_detections(7:end))))

