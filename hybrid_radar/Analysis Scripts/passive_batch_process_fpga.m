function [dec_ref_channel, self_ambg_matrix, cc_matrix] = passive_batch_process_fpga(ref_channel,sur_channel,seg_s,seg_percent,Fs,max_range,exp_dir,zero_padding,td_corr)
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


 %% Decimate segments in to smaller portions
        cc_size = seg_size*(seg_percent/100);  
        dec_ref_channel = seg_ref_channel(1:cc_size,:);
        dec_sur_channel = seg_sur_channel(1:cc_size,:);
        

%% Window Sur channel
        % turn off windowing 
        % for i=1:size(seg_ref_channel,2)
        % dec_sur_channel(:,i) = windowing(dec_sur_channel(:,i), "chebwin");
        % dec_ref_channel(:,i) = windowing(dec_ref_channel(:,i), "chebwin");
        % end
    % Frequency domain implementation
        cc_matrix = complex(zeros((seg_size*zero_padding), (size(ref_channel,1)/seg_size)));
        ref_self_ambg = complex(zeros((seg_size*zero_padding), (size(ref_channel,1)/seg_size)));
    % range limited Xcorr
        for i=1:size(seg_ref_channel,2)
            sur_ch = dec_sur_channel(:,i).';
            ref_ch = dec_ref_channel(:,i).';

            sur_ch=sur_ch/max(real(sur_ch))*16000;
            sur_ch=double(int16(sur_ch));

            ref_ch=ref_ch/max(real(ref_ch))*16000;
            ref_ch=double(int16(ref_ch));
            
            % Compute FFTs
            X1 = fft(sur_ch); %(size(sur_ch,2)*zero_padding)
            X2 = fft(ref_ch);

            X1=X1/128;
            X1=double(int16(X1));

            X2=X2/128;
            X2=double(int16(X2));


            % Compute cross-correlation
            X = X1.*conj(X2);
            % 2,147,483,647
            % 32,767
            X=double(int16((X/2147483647)*32767));
            % figure, hold on,plot(real(X)), plot(imag(X))
            ck = ifft(X,zero_padding*size(X,2));
            ck=double(int16(ck));
            % figure, hold on, plot(real(ck)),plot(imag(ck))
            cc_matrix(:,i) = ck;      

            % make arrays
            r1 = dec_ref_channel(:,i).';

            r1=r1/max(real(r1))*16000;
            r1=double(int16(r1));
            % r2 = dec_ref_channel(:,i).';
            % Compute FFTs
            X3 = fft(r1);
            X3=X3/128;
            X3=double(int16(X3));
            % X4 = fft(r1);
            % Compute cross-correlation
            X5 = X3.*conj(X3);
            X5=double(int16((X5/2147483647)*32767));
            sa = ifft(X5,zero_padding*size(X,2));
            ref_self_ambg(:,i) = sa;    
        end
         [~,bin_zero] = max(ref_self_ambg(:,1));
         cc_matrix = cc_matrix(bin_zero:bin_zero+max_range,:);
         self_ambg_matrix = ref_self_ambg(bin_zero:bin_zero+max_range,:);

      


end


