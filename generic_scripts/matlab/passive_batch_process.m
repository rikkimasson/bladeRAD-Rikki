function cc_matrix = passive_batch_process(ref_channel,sur_channel,seg_s,seg_percent,Fs,max_range)
%PASSIVE_BATCH_PROCESS Summary of this function goes here
%   seg_s : number of segments a second. 
%   seg_percent : percentage of segment used for cross coreclation of 
%                 survallance and reference. Will affect SNR dramatically.


%% Reshape capture into segments
    % Details: segement size determines the limmit of non-ambigious Doppler
    % shift. Comparable to the PRF in active radar. 
    % Non-ambigous Doppler shift = seg_s/2 (Hz)
        seg_size = Fs/seg_s;    % number of samples per segement
        seg_ref_channel = reshape(ref_channel,[seg_size, (size(ref_channel,1)/seg_size)]);
    %seg_sur_channel = reshape(sur_channel,[seg_size, (size(sur_channel,1)/seg_size)]);
        seg_sur_channel = reshape(sur_channel,[seg_size, (size(sur_channel,1)/seg_size)]);
    % plot spectrum of segment of ref and sur channel
        figure
        plot(fftshift(10*log10(abs(fft(seg_ref_channel(:,2000))))))
        hold on 
        plot(fftshift(10*log10(abs(fft(seg_sur_channel(:,2000))))))

 %% Decimate segments in to smaller portions
        cc_size = seg_size*(seg_percent/100);  
        dec_ref_channel = seg_ref_channel(1:cc_size,:);
        dec_sur_channel = seg_sur_channel(1:cc_size,:);
        
%% Window Sur channel
        for i=1:size(seg_ref_channel,2)
        dec_sur_channel(:,i) = windowing(dec_sur_channel(:,i), "Blackman-Harris");
        end

%% Cross-Correlate segments of ref and sur
        cc_matrix = complex(zeros((2*max_range)+1, (size(ref_channel,1)/seg_size)));
    % range limited Xcorr
        for i=1:size(seg_ref_channel,2)
            cc_matrix(:,i) = xcorr(dec_sur_channel(:,i),dec_ref_channel(:,i),max_range); %xcorr(sur_chan,ref_chan) in order to get posative r_bins
        end
        cc_matrix = cc_matrix(max_range+1:end,:); %take zero shifted to +r_max shifted range bins
end

