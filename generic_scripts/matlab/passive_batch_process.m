function [dec_ref_channel, self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,sur_channel,seg_s,seg_percent,Fs,max_range,exp_dir)
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
        freq_axis = linspace(-Fs/2,Fs/2,size(seg_sur_channel,1));
        figure
        fig = plot(freq_axis,fftshift(10*log10(abs(fft(seg_ref_channel(:,1000))))));
        hold on 
        plot(freq_axis,fftshift(10*log10(abs(fft(seg_sur_channel(:,1000))))))
        ylabel('Relative Power (dB)')
        xlabel('Frequency')   
        title("FFT of Single Passive Segment");    
        legend("Reference Channel","Surveillance Channel")
        grid on
        fig_name = exp_dir + "FFT_Segment" + ".jpg";
        saveas(fig,fig_name,'jpeg')

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
        ref_self_ambg = complex(zeros((2*max_range)+1, (size(ref_channel,1)/seg_size)));
    % range limited Xcorr
        for i=1:size(seg_ref_channel,2)
            cc_matrix(:,i) = xcorr(dec_sur_channel(:,i),dec_ref_channel(:,i),max_range);        % xcorr(sur_chan,ref_chan) in order to get positive r_bins
            ref_self_ambg(:,i) = xcorr(dec_ref_channel(:,i),dec_ref_channel(:,i),max_range);   % xcorr
        end
        cc_matrix = cc_matrix(max_range+1:end,:); %take zero shifted to +r_max shifted range bins
        self_ambg_matrix = ref_self_ambg(max_range+1:end,:); %take zero shifted to +r_max shifted range bins
end
