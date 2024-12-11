% load signal and split ref and sur
clear 
clc
close all
addpath('~/Documents/bladeRAD-Rikki/generic_scripts/matlab',...
        '~/Documents/bladeRAD-Rikki/generic_scripts',...
        '~/Documents/bladeRAD-Rikki/generic_scripts/ref_signals/') % path to generic functions



exp_dir= "~/Documents/bladerad_data/Captures/5/"; 
Experiment_ID = 0005;       % Experiment Name



file_location = exp_dir + 'active_' + Experiment_ID;
filename= file_location + ".sc16q11";
% [ref_channel, sur_channel]  = load_passive_data(file_location);
[ signal, signal_i, signal_q ] = load_sc16q11(filename);
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
seg_s = 1000; % number of segments per second - analagos to PRF.
seg_percent = 90;  % percentage of segment used for cross coreclation of 
% survallance and reference. Will affect SNR dramatically.
cc_matrix = passive_batch_process(ref_channel,sur_channel,seg_s,seg_percent,Fs,passive_max_range,exp_dir);
save(exp_dir + 'passive_matrix','cc_matrix')
% RTI Plot
RTI_plot= transpose(10*log10(abs(cc_matrix./max(cc_matrix(:)))));
Range_bin = linspace(0,passive_max_range,size(cc_matrix,1));
time_axis = linspace(0,capture_duration,size(cc_matrix,2));
figure
fig = imagesc(Range_bin,time_axis,RTI_plot,[-50,0]);
% xlim([1 20])
%ylim([0 0.0005])
grid on            
colorbar
ylabel('Time (Sec)')
xlabel('Range Bin')   
fig_title = "Passive RTI - " + Experiment_ID;
title(fig_title);
fig_name = exp_dir + "Passive RTI_" + Experiment_ID + ".jpg";
saveas(fig,fig_name,'jpeg')
saveas(fig,fig_name)

% CAF of entire capture
f_axis = linspace(-seg_s/2,seg_s/2,size(cc_matrix,2));
t_cc_matrix = transpose(cc_matrix);
CAF = fftshift(fft(t_cc_matrix,size(t_cc_matrix,1),1),1);
figure
imagesc(Range_bin,f_axis,10*log10(abs(CAF./max(CAF(:)))),[-100 1]); 
ylim([-500 500])     
% xlim([1 20])
colorbar
ylabel('Doppler Shift (Hz)')
xlabel('Range Bin')  
title("CAF for entire capture" + Experiment_ID)
fig_name = exp_dir + "CAF for entire capture_" + Experiment_ID + ".jpg";
saveas(fig,fig_name,'jpeg')
saveas(fig,fig_name)


% Spectrogram 
        r_bin = 1;
        l_fft = 2056;
        pad_factor = 1;
        overlap_factor = 0.99;
        [spect,f] = spectrogram(cc_matrix(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,seg_s,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-30 0]);   
            ylim([-600 600])
            colorbar
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "Passive Spectrogram - R Bin: " + r_bin + " - " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "Passive Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)



         