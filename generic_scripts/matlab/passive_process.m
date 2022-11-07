function [ref_matrix ,self_ambg_matrix, cc_matrix] = passive_process(ref_channel,sur_channel,xcorr_len,nSurfs,doppler_range,doppler_resolution,Fs,range_bins,exp_dir,range_zero_padding,td_corr)

%PASSIVE_PROCESS Summary of this function goes here

% first determine the truncation factor 
    Ts = 1/Fs;                                % Sample duration
    cap_dur = size(ref_channel,1)/Ts;       % Capture Duration
    xcorr_samps = ceil(xcorr_len/Ts);       % Number of samples to xcorr
    surf_per_sec = cap_dur / nSurfs; % surfaces per second
    cpi_stride = Fs * (1/surf_per_sec);    % number of samples to advance per CPI


% Create doppler shift array and frequency shift matrix
    doppler_bins = doppler_range / doppler_resolution * 2 + 1; % Number of doppler bins
    doppler_shifts = linspace(-doppler_range,doppler_range,doppler_bins); % Doppler shifts for each matched filter

    t = complex((1:xcorr_samps)/Fs); 
    shift_bank = zeros(xcorr_samps,doppler_bins); % bank of refernce segments

    for fD__index = 1:doppler_bins
        shift_bank(:,fD__index) = exp(1j*2*pi*doppler_shifts(fD__index)*t);
    end

% Create cell array for data
    range_doppler_slices = createArrays(nSurfs, [range_bins doppler_bins]);


for i=1:nSurfs

%% truncate CPI in to segments
        seg_ref_channel = ref_channel(1 + ((i-1)*cpi_stride) :xcorr_samps + (i-1)*cpi_stride);
        seg_sur_channel = sur_channel(1 + ((i-1)*cpi_stride) :xcorr_samps + (i-1)*cpi_stride);
    
%% Window Sur channel
        seg_ref_channel = windowing(seg_ref_channel, "Blackman-Harris");
        seg_sur_channel = windowing(seg_sur_channel, "Blackman-Harris");

%% Create bank of matached filters 
        matched_filter_bank = zeros(xcorr_samps,doppler_bins); % bank of refernce segments
      
        for b=1:doppler_bins
        matched_filter_bank(:,b) = seg_ref_channel .* shift_bank(:,b);
        b
        end
%         for  fD__index=1:doppler_bins 
%             fD__index
%             for l=1:xcorr_samps            
%                 matched_filter_bank(:,fD__index) = seg_ref_channel(l)*exp(l*1j*2*pi*doppler_shifts(fD__index)/Fs);
%             end        
%         end

%% Cross-correlate
    range_doppler_slice = zeros(range_bins, doppler_bins);
     parpool(6)
     parfor  fD__index=1:doppler_bins 
            fD__index            
            res = xcorr(seg_sur_channel,matched_filter_bank(:,fD__index),range_bins);
            range_doppler_slice(:,fD__index) = res(range_bins:(2*range_bins)-1);
               
     end

        frame =  10*log(abs(range_doppler_slice./max(range_doppler_slice(:))));

        fig = imagesc(transpose(frame));
        hold on
        ylabel('Doppler Velocity [m/s]')
        xlabel('Range (m)')   
        c = colorbar;
        c.Label.String='Norm Power [dB]';
        colormap jet;
        

     
end


