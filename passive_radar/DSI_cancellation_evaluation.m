
i=85;

p = 0.99;        % subtraction parameter - P must be a positive 
                        % integer less than one to avoid    unwanted discontinuities
                        % arising from zero values in the range–Doppler surface.
thresosthold = 0.001; % cutoff threshold parameter
max_iterations = 500; % maximum number of itterations DSI is CLEANed from CAF slice


% perform CLEAN based DSI Cancellation   
     passive.CLEANed_range_doppler_slices = CLEAN_DSI(passive.range_doppler_slices{i},...
                                                       passive.self_ambg_slices{i},...
                                                       ref_matrix,...
                                                       1,...
                                                       passive.pulses_per_cpi,...
                                                       passive.cpi_stride,...
                                                       passive.zero_padding,...
                                                       passive.range_bins,...
                                                       max_iterations,...
                                                       threshold,p,...
                                                       passive.range_axis,passive.doppler_axis);














% passive.CLEANed_range_doppler_slices


% passive.range_doppler_slices



y_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.doppler_bins);
x_axis = linspace(0,passive.max_range*passive.bin_size*2,passive.range_bins);
x_lim = 500;
y_lim = 100;
dyn_range = inf;

figure
frame =  10*log(abs(passive.range_doppler_slices{i}./max(passive.range_doppler_slices{i}(:))));
% subplot(1,2,1)
        fig = transpose(imagesc(x_axis,y_axis, frame, [-dyn_range 0]));
        hold on
        ylabel('Doppler (Hz)')
        xlabel('Bistatic Range (m)')   
        c = colorbar;
        c.Label.String='Norm Power (dB)';
        ylim([-y_lim y_lim]);
        xlim([-inf x_lim]);
        colormap jet;
subplot(1,2,2)  
figure
frame =  10*log(abs(passive.CLEANed_range_doppler_slices{i}./max(passive.CLEANed_range_doppler_slices{i}(:))));

        fig = transpose(imagesc(x_axis,y_axis, frame, [-dyn_range 0]));
        hold on
        ylabel('Doppler (Hz)')
        xlabel('Bistatic Range (m)')   
        c = colorbar;
        c.Label.String='Norm Power (dB)';
        ylim([-y_lim y_lim]);
        xlim([-inf x_lim]);
        colormap jet;
