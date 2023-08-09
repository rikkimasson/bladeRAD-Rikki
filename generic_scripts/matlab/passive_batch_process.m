function [dec_ref_channel, self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,sur_channel,seg_s,seg_percent,Fs,max_range,exp_dir,zero_padding,td_corr)
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
%         fig_name = exp_dir + "FFT_Segment" + ".jpg";
%         saveas(fig,fig_name,'jpeg')

 %% Decimate segments in to smaller portions
        cc_size = seg_size*(seg_percent/100);  
        dec_ref_channel = seg_ref_channel(1:cc_size,:);
        dec_sur_channel = seg_sur_channel(1:cc_size,:);
        
%% Cross-Correlate segments of ref and sur
if td_corr == true
'Time Domain Based XCORR'
% Time domain implementation
        cc_matrix = complex(zeros((zero_padding*(2*max_range+1)), (size(ref_channel,1)/seg_size)));
        ref_self_ambg = complex(zeros((zero_padding*(2*max_range+1)), (size(ref_channel,1)/seg_size)));
    % range limited Xcorr
        for i=1:size(seg_ref_channel,2)
            if zero_padding > 1
%                 int_ref_channel = ifft([ fft(dec_ref_channel(:,i)); zeros(size(dec_ref_channel,1)*(zero_padding-1),1) ]);
%                 int_sur_channel = ifft([ fft(dec_sur_channel(:,i)); zeros(size(dec_ref_channel,1)*(zero_padding-1),1) ]);
                cor_prod = xcorr(dec_sur_channel(:,i),dec_ref_channel(:,i),max_range);
                cc_matrix(:,i) = ifft( [ fft(cor_prod); zeros(size(cor_prod,1)*(zero_padding-1),1) ]);             % xcorr(sur_chan,ref_chan) in order to get positive r_bins
                cor_prod = xcorr(dec_ref_channel(:,i),dec_ref_channel(:,i),max_range);
                ref_self_ambg(:,i) = ifft( [ fft(cor_prod); zeros(size(cor_prod,1)*(zero_padding-1),1) ]);

            else
            cc_matrix(:,i) = xcorr(dec_sur_channel(:,i),dec_ref_channel(:,i),max_range);       % xcorr(sur_chan,ref_chan) in order to get positive r_bins
            ref_self_ambg(:,i) = xcorr(dec_ref_channel(:,i),dec_ref_channel(:,i),max_range);   % xcorr
            end
        end
         [~,bin_zero] = max(cc_matrix(:,1));
         cc_matrix = cc_matrix(bin_zero:bin_zero+max_range,:);
         self_ambg_matrix = ref_self_ambg(bin_zero:bin_zero+max_range,:);
         return
else
'Frequency Based XCORR'

%% Window Sur channel
        for i=1:size(seg_ref_channel,2)
        dec_sur_channel(:,i) = windowing(dec_sur_channel(:,i), "chebwin");
        dec_ref_channel(:,i) = windowing(dec_ref_channel(:,i), "chebwin");
        end
    % Frequency domain implementation
        cc_matrix = complex(zeros((seg_size*zero_padding), (size(ref_channel,1)/seg_size)));
        ref_self_ambg = complex(zeros((seg_size*zero_padding), (size(ref_channel,1)/seg_size)));
    % range limited Xcorr
        for i=1:size(seg_ref_channel,2)
            
            % make arrays
                %sur_ch = [ complex(zeros(seg_size-1,1)).' dec_sur_channel(:,i).' ];
                %ref_ch = [ dec_ref_channel(:,i).' complex(zeros(seg_size-1,1)).'  ];
%                 sur_ch = [ dec_sur_channel(:,i).' complex(zeros((seg_size*zero_padding),1)).' ];
%                 ref_ch = [ complex(zeros((seg_size*zero_padding),1)).'  dec_ref_channel(:,i).'  ];
%                 sur_ch = [ dec_sur_channel(:,i).' complex(zeros((seg_size),1)).' ];
%                 ref_ch = [ complex(zeros((seg_size),1)).'  dec_ref_channel(:,i).'  ];
                    sur_ch = dec_sur_channel(:,i).';
                    ref_ch = dec_ref_channel(:,i).';


            % Compute FFTs
            X1 = fft(sur_ch); %(size(sur_ch,2)*zero_padding)
            X2 = fft(ref_ch);
            % Compute cross-correlation
            X = X1.*conj(X2);
            ck = ifft(X,zero_padding*size(X,2));
            cc_matrix(:,i) = ck;      

            % make arrays
            r1 = dec_ref_channel(:,i).';
            r2 = dec_ref_channel(:,i).';
            % Compute FFTs
            X3 = fft(r1);
            X4 = fft(r2);
            % Compute cross-correlation
            X5 = X3.*conj(X4);
            sa = ifft(X5,zero_padding*size(X,2));
            ref_self_ambg(:,i) = sa;    
        end
         [~,bin_zero] = max(ref_self_ambg(:,1));
         cc_matrix = cc_matrix(bin_zero:bin_zero+max_range,:);
         self_ambg_matrix = ref_self_ambg(bin_zero:bin_zero+max_range,:);

      


end


