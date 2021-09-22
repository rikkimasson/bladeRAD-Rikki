%% Capture Parameters
C = physconst('LightSpeed');
Bw = 30e6;         
Fs=2*Bw;
Cap_dur = 5
Samp_dur = 1 / Bw
Num_cap_samps = Cap_dur * 2 / Samp_dur
f_size = Num_cap_samps * 16 

%% Radar Parameters 1
Fc = 5.22e9;
Rx1_gain = 50;
Rx2_gain = 50;

Test_id = 1211;

RF_freq = Fc/1e6; %in MHz 
Bw_M = Bw/1e6;

save_directory = "/media/piers/PassiveRadar/15-04-21/" + Test_id;
sys_command = "mkdir " + save_directory;
system(sys_command)

%SDR 2
command = "/home/piers/Documents/bladeRF_Code/PassiveRadar/Shell_Scripts/PassiveRadar.sh " + Num_cap_samps + " " + Rx1_gain + " " + " " + Rx2_gain + " " + RF_freq + " " + Bw_M + " " + Test_id; 

status = system(command);


r_max = 10;

%% Load Signal    
    %load rx file
        tic
%         rx_file= "/media/piers/data_drive/BladeRF_Experiments/Captures/Passive/" + Test_id + "/" + Test_id + ".sc16q11"
%         rx_file = "/Volumes/data_drive/BladeRF_Experiments/Captures/Passive/" + Test_id + "/" + Test_id + ".sc16q11"
        rx_file = "/media/piers/PassiveRadar/14-04-21/" + Test_id + "/" + Test_id + ".sc16q11"
%         rx_file= Test_id + "\" + Test_id + ".sc16q11"
        
        raw_data= load_sc16q11(rx_file); 
        loading_rx_file = toc 
    
    % De-interleave data 
        ref_channel  = raw_data(1:2:end,:);
        sur_channel  = raw_data(2:2:end,:);
        clear raw_data
%         save 5701.mat ref_channel sur_channel Fs Test_id -v7.3
    % Plot time domain signals
         figure
         plot(real(ref_channel))
         title("Ref channel time series");
         figure
         plot(real(sur_channel))
         title("Sur channel time series");




%% Psuedo Wi-Fi signal
        %     ps_sur_channel = ref_channel;
        %     ps_sur_channel = circshift(ps_sur_channel,10); % shift by 10 range bins
        %     for i=1:Num_cap_samps/2
        %             ps_sur_channel(i) = ps_sur_channel(i) * exp(i * 1j * 2 * pi * 200 / Fs);
        %             % Range loop 
        %     end

%% Reshape capture into segments
    %Details: segement size determines the limmit of non-ambigious Doppler
    %shift. Comparable to the PRF in active radar. 
    %Non-ambigous Doppler shift = seg_s/2 (Hz)
        seg_s = 5000;           % number of segments a second. 
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
        seg_percent = 20;           % percentage of segment used for cross coreclation of survallance and reference
        cc_size = seg_size*10/100;  
        dec_ref_channel = seg_ref_channel(1:cc_size,:);
        dec_sur_channel = seg_sur_channel(1:cc_size,:);

 %% DSI Cancellation      
%      % The required time domain filter parameters for each filter are passed 
%      % through the kwargs interface class
%         tdfp = kwargs;
%         tdfp.K = 10;          % multipath echoes are backscattered from the first K range bins
%         tdfp.D = 512;     % Maximum Doppler extension measured in Doppler bins
% %         tdfp.T      % Number of batches for ECA-B or batch length for block LMS
% %         tdfp.Na     % Sliding window size, measured in samples
% %         tdfp.imp    % Implementation type
% %         tdfp.mu     % Step size parameter for the iterative algorithms - LMS or NLMS
% %         tdfp.lamb   % Forgetting factor for the RLS algoritm
% %         tdfp.ui     % Update interval for the iterative algorithms - LMS, NLMS, RLS
% %         tdfp.w_init % Initialization vector for the iterative algorithms - LMS, NLMS, RLS, (default: None)
% 
%         filtered_surv_ch = complex(zeros(size(dec_ref_channel)));
%         for i=1:size(seg_ref_channel,2)
%         filtered_surv_ch(:,i) = time_domain_filter_surveillance(dec_ref_channel(:,i), dec_sur_channel(:,i), "ECA", tdfp);
%         end
%         toc       
%  
%         
        
%% Window Sur channel
        for i=1:size(seg_ref_channel,2)
        dec_sur_channel(:,i) = windowing(dec_sur_channel(:,i), "Blackman-Harris");
        end
        toc
%% Cross-Correlate segments of ref and sur
        cc_matrix = complex(zeros((2*r_max)+1, (size(ref_channel,1)/seg_size)));
    % range limited Xcorr
        tic
        for i=1:size(seg_ref_channel,2)
        cc_matrix(:,i) = xcorr(dec_sur_channel(:,i),dec_ref_channel(:,i),r_max); %xcorr(sur_chan,ref_chan) in order to get posative r_bins
        end
        toc
        cc_matrix = cc_matrix(r_max+1:end,:); %take zero shifted to +r_max shifted range bins


%% RTI Plot
    RTI_plot= transpose(10*log10(abs(cc_matrix./max(cc_matrix(:)))));

    Range_bin = linspace(0,r_max,size(cc_matrix,1));
    time_axis = linspace(0,Cap_dur,size(cc_matrix,2));

    figure
    fig = imagesc(Range_bin,time_axis,RTI_plot,[-50,0]);
    % xlim([1 20])
    %ylim([0 0.0005])
    grid on            
    colorbar
    ylabel('Time (Sec)')
    xlabel('Range Bin')   
    fig_title = "Psudo Range Shifted";
    title(fig_title);
    fig_name = save_directory + "/RTI_" + Test_id + ".jpg";
    saveas(fig,fig_name,'jpeg')
    plot_signal = toc
    % 

%% CAF of entire capture
    f_axis = linspace(-seg_s/2,seg_s/2,size(cc_matrix,2));
    t_cc_matrix = transpose(cc_matrix);
    CAF = fftshift(fft(t_cc_matrix,size(t_cc_matrix,1),1),1);
    figure
    imagesc(Range_bin,f_axis,10*log10(abs(CAF./max(CAF(:)))),[-50 1]); 
    ylim([-500 500])     
    % xlim([1 20])
    colorbar
     ylabel('Doppler Shift (Hz)')
     xlabel('Range Bin')   
     fig_title = "Psudo Range % Doppler Shifted";
 

%% Spectrogram 
    r_bin = 1;
    l_fft = 1024;
    pad_factor = 1;
    overlap_factor = 0.99;
    [spect,f] = spectrogram(cc_matrix(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,seg_s,'centered','yaxis');
    % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
    v=dop2speed(f,C/Fc)*2.237;
    spect= 10*log10(abs(spect./max(spect(:))));
    figure
    fig = imagesc(time_axis,f,spect,[-50 0]);   
    ylim([-600 600])
    colorbar
    xlabel('Time (Sec)')
    % ylabel('Radial Velocity (mph)')   
    ylabel('Doppler Frequency (Hz)')  
    fig_title = "Spectrogram - R Bin: " + r_bin + " - Test id: " + Test_id;
    title(fig_title);
    fig_name = save_directory + "/Spectrogram_" + Test_id + ".jpg";
    saveas(fig,fig_name,'jpeg')
    % 

% %% MTI Filtering 
% % Single Delay Line Filter 
% MTI_Data = zeros(size(cc_matrix));
%       for i=2:size(cc_matrix,2)
%             MTI_Data(:,i) = cc_matrix(:,i)-cc_matrix(:,i-1);
%       end
%       
% % %Plot MTI RTI      
% % MTI_RTI_plot= transpose(10*log10(abs(MTI_Data./max(MTI_Data(:)))));
% % figure
% % fig = imagesc(Range_bin,time_axis,MTI_RTI_plot,[-50,0]);
% % xlim([1 20])
% % %ylim([0 0.0005])
% % grid on            
% % colorbar
% % ylabel('Time (Sec)')
% % xlabel('Range Bin')   
% % fig_title = "Monostatic Single Delay Line MTI  RTI - Test " + Test_id;
% % title(fig_title);
% % fig_name = save_directory + "/MTI_RTI_" + Test_id + ".jpg";
% % saveas(fig,fig_name,'jpeg')
% % plot_signal = toc     
% 
% % 
% % %Plot MTI Spectrogram  
% [spect,f] = spectrogram(MTI_Data(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,seg_s,'centered','yaxis');
% % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
% v=dop2speed(f,C/Fc)*2.237;
% spect= 10*log10(abs(spect./max(spect(:))));
% figure
% fig = imagesc(time_axis,f,spect,[-30 0]);
% ylim([-600 600])
% colorbar
% xlabel('Time (Sec)')
% % ylabel('Radial Velocity (mph)')   
% ylabel('Doppler Frequency (Hz)')  
% fig_title = "Monostatic Single Delay Line MTI Spectrogram - Test " + Test_id;
% title(fig_title);
% fig_name = save_directory + "/MTI_Spectrogram_" + Test_id + ".jpg";
% saveas(fig,fig_name,'jpeg')
% 
