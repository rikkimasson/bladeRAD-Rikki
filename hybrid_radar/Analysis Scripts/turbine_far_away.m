tic
% 1. Respond in the command window when prompted in order to point MATLAB
% to the correct directory
% prompt = "Which computer is running the script? 1=Rikki's Mac, 2=Jacob's Laptop, 3=Rikki UCL \n"
clear
clc
close all
pc = 3;%input(prompt);
addpath("C:/Users/rw/Documents/RADCOM/misc_helper_functions");

local_save_directory = "D:/passive_archery_data/archery_wind_turbine_32x/";
repo_directory = "C:\Users\rw\Documents\";


% 2. Processing Flags - choose whether to
process_active_flag = false;
process_passive_flag = true;
ground_truth_flag = false;

% 3. Select the experiment you wish to process


% adds paths to generic functions

addpath(repo_directory + '\bladeRAD-Rikki\generic_scripts\matlab',...
repo_directory + '\bladeRAD-Rikki\generic_scripts\matlab\CFAR\',...
repo_directory + '\bladeRAD-Rikki\generic_scripts',...
repo_directory + '\bladeRAD-Rikki\generic_scripts\ref_signals\')
  


%% Load .mat file containing experiment setup parameters

    % mat_file_name = local_save_directory + i + bracket + "Experimental Configuration.mat";
    % load(mat_file_name);
    exp_dir = local_save_directory ;


%% Passive Processing

    if process_passive_flag == true
        passive.min_range = 1; % maximum number of range bins - xcorr shifts.
        passive.max_range = 500; % maximum number of range bins - xcorr shifts.
        passive.range_zero_padding = 1; % 1 = none, 2 = 100%
        passive.td_corr = false; % true = time domain xcorr; false = freq domain xcorr
        passive.seg_s = 1000; % number of segments per second - analagos to PRF.
        passive.seg_percent = 100; % % of segment used for xcorr.
        
        % load signal and split ref and sur
        file_location = local_save_directory + 'passive_data_test_1/';
        load(file_location+'passive_ref_output');
        data=double(data);
        ref_channel  = data(1:2:end)+1j .*data(2:2:end);

        file_location = local_save_directory + 'passive_data_test_1/';
        load(file_location+'passive_surv_output');
        data=double(data);
        
        sur_channel  = data(1:2:end)+1j .*data(2:2:end);
        clearvars 'data' 
       
        
        % Add Noise to Passive Radar
        if (0)
        noisy_sur_channel= awgn(sur_channel,-2,'measured');
        else
            noisy_sur_channel=sur_channel;
        end

        figure
        plot(abs(fft(ref_channel(1:100000))))

        myx=xcorr(ref_channel(1:400000),ref_channel(1:400000));

        figure,plot(abs(myx));

        % myxcorr(ref_channel(1:400000),ref_channel(1:400000));

        % mresult=ref_channel(1:400000).*conj(ref_channel(1+55440:400000+55440));

        % figure,
        % hold on
        % plot(real(mresult))
        % plot(imag(mresult))
        % plot(abs(mresult))
        % % 
        % figure
        % for h=1:1:2000
        %     mstart=h;
        %     mspacing=8192;
        %     XF=fft(ref_channel(mstart:mstart+mspacing));
        %     scatter(real(XF),imag(XF))
        % 
        % end

        % figure
        % plot(abs(fft(ref_channel(1:10000))))
        % 
        %  figure
        % plot(real((ref_channel(:))))
        % 
        % figure,plot(abs(xcorr(ref_channel(1:20000),ref_channel(1:20000))))
        % 
        % 
        % figure,plot(abs(xcorr(ref_channel(1:20000),ref_channel(1:20000),'unbiased')))
        % mshift=2986;
        % mout=sum(abs((ref_channel(1:20000).*conj(ref_channel(1+mshift:20000+mshift)))));
        % 
        % [out] = myxcorr(ref_channel(1:20000),ref_channel(1:20000));
        % 
        % 
        % figure,plot(abs(mout))
        % 
        % figure
        % hold on
        % plot(real(ref_channel(1:20000)))
        % plot(real((ref_channel(1+2987:20000+2987))))
        % 
        % 
        % temp3=zeros(1,40000);
        % for i=1:1:40000
        %     temp3(i)=(sum(ref_channel(i:i+400).*conj(ref_channel(i+2987:i+400+2987))));    
        % end
        % figure,plot(abs(temp3))
        % 
        % figure,plot(angle(temp3))
        % 
        % XF=fft(ref_channel(17994:17994+2987));
        % 
        % figure
        % scatter(real(XF(1:10)),imag(XF(1:10)))

        % XX=ifft(fft(ref_channel).*conj(fft(noisy_sur_channel)));
        % 
        % XX_auto=ifft(fft(ref_channel).*conj(fft(ref_channel)));
        % 
        % figure
        % plot(20*log(abs(XX(1:10000))))
        % 
        % figure
        % hold on
        % plot(20*log(abs(XX(1:10000))))
        % plot(20*log(abs(XX_auto(1:10000))))
        % 


        % [ref_matrix_fpga ,self_ambg_matrix_fpga, cc_matrix_fpga] = passive_batch_process_fpga(ref_channel,noisy_sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.max_range,exp_dir,passive.range_zero_padding,passive.td_corr);
       

        % self_ambg_matrix=self_ambg_matrix_fpga;
        % cc_matrix=cc_matrix_fpga;
        % ref_matrix= ref_matrix_fpga;
        passive.Fs=15e6;
        C=299792458;
        passive.Fc=498e6;
        capture_duration=(length(ref_channel)/passive.Fs);%15;
        ref_channel=ref_channel';
        noisy_sur_channel=noisy_sur_channel';

        [ref_matrix ,self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,noisy_sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.min_range,passive.max_range,exp_dir,passive.range_zero_padding,passive.td_corr);
        % save(exp_dir + 'passive_matrix','cc_matrix')


        
%% Proccess Passive data into Range-Doppler Slices

        passive.PRF = passive.seg_s; % seg_s
        passive.cpi = 2.5; % coherent proccessing interval (s)
        passive.cpi_overlap = 0.8; % overlap between CPIs (watch this - too large will cause slow exceution)
        passive.doppler_window = 'Blackman-Harris';
        passive.dopp_zero_padding = 1;
        passive.dynamic_range = +50; %dynamic range of image
        passive.range_bin_size = (1/passive.Fs)/passive.range_zero_padding * C/2;
        passive.max_range_m = passive.max_range*passive.range_bin_size;
        % Create range-Doppler surfaces
        [passive.number_cpi,...
            passive.pulses_per_cpi,...
            passive.range_doppler_slices] = rangeDopplerSlice(cc_matrix,passive.cpi,passive.PRF,...
            passive.cpi_overlap,...
            passive.dopp_zero_padding,...
            passive.doppler_window);
        % Create Ambiguity surfaces
        [~,~,passive.self_ambg_slices] = rangeDopplerSlice(self_ambg_matrix,passive.cpi,passive.PRF,...
            passive.cpi_overlap,...
            passive.dopp_zero_padding,...
            passive.doppler_window);
        
        % Create Active Data Range and Doppler Axis
        
        passive.cpi_stride = round(passive.pulses_per_cpi*(1-passive.cpi_overlap)); % number of pulses to stride each for next CPI
        passive.velocity_conv = C*(((1/C)/(passive.Fc/C)));
        passive.no_range_bins = size(cc_matrix,1);
        passive.range_bins = 1:passive.no_range_bins;
        passive.doppler_bins = passive.pulses_per_cpi*passive.dopp_zero_padding+1;
        passive.doppler_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.doppler_bins);
        passive.doppler_velocity_axis = passive.doppler_axis*passive.velocity_conv;
        passive.range_axis = linspace(0,passive.max_range_m,passive.no_range_bins);
        
        
        % create video of Passive range-Doppler slices        
        video_name = exp_dir + "passive_range-Doppler_" + 'turbine' + ".avi";
        video_title = "Passive Pre-DSI";
        dynamic_range = +inf;
        max_range = 11000;
        max_doppler = 50;
        frame_rate = 1/(capture_duration/passive.number_cpi);
        createVideo(passive.range_doppler_slices,frame_rate,...
            passive.range_axis,max_range,...
            passive.doppler_axis,max_doppler,...
            dynamic_range,video_name,video_title);
        %
%% Direct Signal Interference Cancellation
% set DSI cancellation parameters

        p = 0.999;        % subtraction parameter - P must be a positive
        % integer less than one to avoid unwanted discontinuities
        % arising from zero values in the range–Doppler surface.
        threshold = 0.005; % cutoff threshold parameter
        max_iterations = 100; % maximum number of itterations DSI is CLEANed from CAF slice
        number_rbins = size(passive.range_doppler_slices{1},2);
        
        % perform CLEAN based DSI Cancellation
        passive.CLEANed_range_doppler_slices = CLEAN_DSI(passive.range_doppler_slices,...
            passive.self_ambg_slices,...
            ref_matrix,...
            passive.number_cpi,...
            passive.pulses_per_cpi,...
            passive.cpi_stride,...
            passive.dopp_zero_padding,...
            number_rbins,...
            max_iterations,...
            threshold,p,...
            passive.range_axis,passive.doppler_axis);
        
        %   % create video of CLEANed range-Doppler slices
        video_name = exp_dir + "CLEANed_range-Doppler_CLEANed_log_Exp_" + '_turbine' + ".avi";
        %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";
        video_title = "CLEANed Passive Radar Capture";
        dynamic_range = 50;
        max_range = 5000;
        max_doppler = 50;
        frame_rate = 1/(capture_duration/passive.number_cpi);
        createVideo(passive.CLEANed_range_doppler_slices,frame_rate,...
            passive.range_axis,max_range,...
            passive.doppler_velocity_axis,max_doppler,...
            dynamic_range,video_name,video_title);        
    end 

toc


% temp_full=passive.CLEANed_range_doppler_slices{45};

% save('image_frame_save','temp_full')

temp=passive.CLEANed_range_doppler_slices{45};

load('image_frame_save')

figure
frame = 20*log10(abs(temp./max(temp(:))));
fig = transpose(imagesc(passive.range_axis,passive.doppler_velocity_axis,frame, [-dynamic_range 0]));
clim([-dynamic_range 0]) % -dynamic_range
c = colorbar;

figure
frame_og = 20*log10(abs(temp_full./max(temp_full(:))));
fig = transpose(imagesc(passive.range_axis,passive.doppler_velocity_axis,frame_og, [-dynamic_range 0]));
clim([-dynamic_range 0]) % -dynamic_range
c = colorbar;

diff_frame = 20*log10(abs((temp_full-temp)./max(temp_full(:))));
figure
fig = transpose(imagesc(passive.range_axis,passive.doppler_velocity_axis,diff_frame, [-dynamic_range 0]));
clim([-dynamic_range 0]) % -dynamic_range
c = colorbar;



