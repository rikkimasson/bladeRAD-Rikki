%% Fuse active and passive range Doppler images 

%% Interpolate passive radar data
    range_factor = passive.range_bin_size/active.range_bin_size; % Interolation factor in range
    doppler_factor = passive.velocity_conv/active.velocity_conv; % Interpolation factor in doppler 
    
    passive.inter_range_doppler_slices = interpolateRangeDoppler(passive.CLEANed_range_doppler_slices,range_factor,doppler_factor);
  
    [passive.interp_doppler_bins passive.interp_range_bins] = size(passive.inter_range_doppler_slices{1});
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
     createVideo(passive.inter_range_doppler_slices,frame_rate,...
                 passive.range_axis,max_range,...
                 passive.doppler_axis,max_doppler,...
                 dynamic_range,video_name,video_title);


%% Scale range Doppler images and fuse
passive_range_slide = 0;
active_range_slide  = 2; %2
[nAdbins nArbins] = size(active.range_doppler_slices{1});
[nPdbins nPrbins] = size(passive.inter_range_doppler_slices{1});
fused_range = min([nArbins nPrbins])-max([active_range_slide passive_range_slide]);
fused_doppler = (min([nAdbins nPdbins]))/2-1;
hybrid.range_doppler_slices = createArrays(size(active.range_doppler_slices,2), [fused_range ceil(fused_doppler)]);

wv = 'db2';
lv = 5;
active_data = createArrays(size(active.range_doppler_slices,2), [fused_range ceil(fused_doppler)]);
passive_data = createArrays(size(active.range_doppler_slices,2), [fused_range ceil(fused_doppler)]);

for i=1:size(active.range_doppler_slices,2)
    active_data{i} = flipud(active.range_doppler_slices{i}(nAdbins/2-fused_doppler:nAdbins/2+fused_doppler,1+active_range_slide:fused_range+active_range_slide));
    passive_data{i} = passive.inter_range_doppler_slices{i}(nPdbins/2-fused_doppler:nPdbins/2+fused_doppler,1+passive_range_slide:fused_range+passive_range_slide);
    active_data1 = (active_data{i}./max(active_data{i}(:)));
    passive_data1 = (passive_data{i}./max(passive_data{i}(:)));
%     hybrid.range_doppler_slices{i} = wfusmat(abs(active_data{i}),abs(passive_data{i}),'mean');
    [aA,aH,aV,aD] = dwt2(abs(active_data{i}),'sym5','mode','per');
    [pA,pH,pV,pD] = dwt2(abs(passive_data{i}),'sym5','mode','per');
    
    hybrid.range_doppler_slices{i} = idwt2((aA),(aH+pH),(aV+pV),(aD+pD),'db4')

%     hybrid.range_doppler_slices{i} = abs(active_data1) + abs(passive_data1);

end 

hybrid.range_axis = active.range_axis(1:fused_range);
hybrid.doppler_velocity_axis = active.doppler_velocity_axis(nAdbins/2-fused_doppler:nAdbins/2+fused_doppler);

% % create video of  range-Doppler slices
% %      video_name = exp_dir + "Fused_Range_Doppler_Slices" + Experiment_ID + ".avi";
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
    dynamic_range = +90;
    max_range = +inf;
    max_doppler = 60;
    createVideo3(active_data, ...
                passive_data,hybrid.range_doppler_slices, ...
                frame_rate,hybrid.range_axis,hybrid.doppler_velocity_axis, ...
                max_range,max_doppler,dynamic_range,video_name,video_title)



%% SNR of peak return 
passive_range_lookup = ones(size(passive.CLEANed_range_doppler_slices{2})).*passive.range_axis;
passive_doppler_lookup = ones(size(passive.CLEANed_range_doppler_slices{2})).*transpose(passive.doppler_velocity_axis);

passive_range_lookup = ones(size(passive.inter_range_doppler_slices {2})).*passive.interp_range_axis;
passive_doppler_lookup = ones(size(passive.inter_range_doppler_slices{2})).*transpose(passive.interp_doppler_velocity_axis);

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

   [target_mag, index] = max(passive.inter_range_doppler_slices{i}(:));
   passive_range_detections(i) = passive_range_lookup(index);
   passive_target_power = 10*log10(abs(target_mag));
   passive_noise_estimate = 10*log10(mean(abs(passive.inter_range_doppler_slices{i}(350:500,110:120)),'all'));   
   passive_detection_snr(i) = passive_target_power - passive_noise_estimate;


   [target_mag, index] = max(active.range_doppler_slices{i}(:));
   active_range_detections(i) = active_range_lookup(index);
   active_target_power = 10*log10(abs(target_mag));
   active_noise_estimate = 10*log10(mean(abs(active.range_doppler_slices{i}(350:500,110:120)),'all'));
   active_detection_snr(i) = active_target_power - active_noise_estimate;

   [target_mag, index] = max(hybrid.range_doppler_slices{i}(:));
   hybrid_range_detections(i) = hybrid_range_lookup(index);
   hybrid_target_power = 10*log10(abs(target_mag));
   hybrid_noise_estimate = 10*log10(mean(abs(hybrid.range_doppler_slices{i}(:,110:120)),'all'));
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
