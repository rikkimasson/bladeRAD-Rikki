clear
clc
close all
tic
% 1. Respond in the command window when prompted in order to point MATLAB
% to the correct directory
% prompt = "Which computer is running the script? 1=Rikki's Mac, 2=Jacob's Laptop, 3=Rikki UCL \n"
pc = 3;%input(prompt);
addpath("C:/Users/rw/Documents/RADCOM/misc_helper_functions");

local_save_directory = "D:/passive_archery_data/run_5_32x_small/";
local_save_directory = "C:/Users/rw/Downloads/run_5_32x_small/run_5_32x_small/";
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
        passive.min_range = 2; % maximum number of range bins - xcorr shifts.
        passive.max_range = 20; % maximum number of range bins - xcorr shifts.
        passive.range_zero_padding = 1; % 1 = none, 2 = 100%
        passive.td_corr = false; % true = time domain xcorr; false = freq domain xcorr
        passive.seg_s = 2000; % number of segments per second - analagos to PRF.
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

        % write_binary_to_text('passive_ref',ref_channel(1:2000000),16,0);
        % 
        % write_binary_to_text('passive_surv',sur_channel(1:2000000),16,0);
    %%
permanant_carriers= [0 48 54 87 141 156 192 ...
201 255 279 282 333 432 450 ...
483 525 531 618 636 714 759 ...
765 780 804 873 888 918 939 ...
942 969 984 1050 1101 1107 1110 ...
1137 1140 1146 1206 1269 1323 1377 ...
1491 1683 1704 1752 1758 1791 1845 ...
1860 1896 1905 1959 1983 1986 2037 ...
2136 2154 2187 2229 2235 2322 2340 ...
2418 2463 2469 2484 2508 2577 2592 ...
2622 2643 2646 2673 2688 2754 2805 ...
2811 2814 2841 2844 2850 2910 2973 ...
3027 3081 3195 3387 3408 3456 3462 ...
3495 3549 3564 3600 3609 3663 3687 ...
3690 3741 3840 3858 3891 3933 3939 ...
4026 4044 4122 4167 4173 4188 4212 ...
4281 4296 4326 4347 4350 4377 4392 ...
4458 4509 4515 4518 4545 4548 4554 ...
4614 4677 4731 4785 4899 5091 5112 ...
5160 5166 5199 5253 5268 5304 5313 ...
5367 5391 5394 5445 5544 5562 5595 ...
5637 5643 5730 5748 5826 5871 5877 ...
5892 5916 5985 6000 6030 6051 6054 ...
6081 6096 6162 6213 6219 6222 6249 ...
6252 6258 6318 6381 6435 6489 6603 ...
6795 6816];

permanant_carriers=permanant_carriers+1;

%% mapping from 64qam to bits
map_64qam=cell(8,8);
map_64qam{1,1}='100000';
map_64qam{1,2}='100010';
map_64qam{1,3}='101010';
map_64qam{1,4}='101000';
map_64qam{1,5}='001000';
map_64qam{1,6}='001010';
map_64qam{1,7}='000010';
map_64qam{1,8}='000000';

map_64qam{2,1}='100001';
map_64qam{2,2}='100011';
map_64qam{2,3}='101011';
map_64qam{2,4}='101001';
map_64qam{2,5}='001001';
map_64qam{2,6}='001011';
map_64qam{2,7}='000011';
map_64qam{2,8}='000001';

map_64qam{3,1}='100101';
map_64qam{3,2}='100111';
map_64qam{3,3}='101111';
map_64qam{3,4}='101101';
map_64qam{3,5}='001101';
map_64qam{3,6}='001111';
map_64qam{3,7}='000111';
map_64qam{3,8}='000101';

map_64qam{4,1}='100100';
map_64qam{4,2}='100110';
map_64qam{4,3}='101110';
map_64qam{4,4}='101100';
map_64qam{4,5}='001100';
map_64qam{4,6}='001110';
map_64qam{4,7}='000110';
map_64qam{4,8}='000100';

map_64qam{5,1}='110100';
map_64qam{5,2}='110110';
map_64qam{5,3}='111110';
map_64qam{5,4}='111100';
map_64qam{5,5}='011100';
map_64qam{5,6}='011110';
map_64qam{5,7}='010110';
map_64qam{5,8}='010100';

map_64qam{6,1}='110101';
map_64qam{6,2}='110111';
map_64qam{6,3}='111111';
map_64qam{6,4}='111101';
map_64qam{6,5}='011101';
map_64qam{6,6}='011111';
map_64qam{6,7}='010111';
map_64qam{6,8}='010101';

map_64qam{7,1}='110001';
map_64qam{7,2}='110011';
map_64qam{7,3}='111011';
map_64qam{7,4}='111001';
map_64qam{7,5}='011001';
map_64qam{7,6}='011011';
map_64qam{7,7}='010011';
map_64qam{7,8}='010001';

map_64qam{8,1}='110000';
map_64qam{8,2}='110010';
map_64qam{8,3}='111010';
map_64qam{8,4}='111000';
map_64qam{8,5}='011000';
map_64qam{8,6}='011010';
map_64qam{8,7}='010010';
map_64qam{8,8}='010000';

%%      
        
% Add Noise to Passive Radar
if (0)
noisy_sur_channel= awgn(sur_channel,-2,'measured');
else
    noisy_sur_channel=sur_channel;
end

figure,plot(real(data(1:100000)))

figure,plot(abs(fft(ref_channel(1:100000))))

fs=3.84e9/8/32;
dt=1/fs;
delta=[224, 112,56,28]*1e-6;
dd=delta(4);
Tu=896e-6;
D=dd/dt;
D=413;
S=Tu/dt;
Tots=S+D;
B=S-D;
P1=zeros(1,100000);

temp=ref_channel(1:200000);

%% try this again
numb_carriers=6817;
fs=3.84e9/8/32;
dt=1/fs;
delta=[224, 112,56,28]*1e-6;
dd=delta(4);
Tu=896e-6;
D=dd/dt;
D=413;
S=Tu/dt;
Tots=S+D;
B=S-D;
% [outputArg1] = produce_ideal_ofdm_symbol(ref_channel(1:200000),S,D,numb_carriers,permanant_carriers,dt);
[ref_channel_ideal] = produce_ideal_ofdm_symbol(ref_channel,S,D,numb_carriers,permanant_carriers,dt);

% figure,
% hold on
% plot(real(outputArg1))
% plot(imag(outputArg1))
%% back to regularly schedualed programming 
% 
passive.Fs=15e6;
C=299792458;
passive.Fc=498e6;
capture_duration=(length(ref_channel)/passive.Fs);%15;
% ref_channel=ref_channel';
ref_channel=ref_channel_ideal';
noisy_sur_channel=noisy_sur_channel';

[ref_matrix ,self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,noisy_sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.min_range,passive.max_range,exp_dir,passive.range_zero_padding,passive.td_corr);
% save(exp_dir + 'passive_matrix','cc_matrix')


        
%% Proccess Passive data into Range-Doppler Slices

        passive.PRF = passive.seg_s; % seg_s
        passive.cpi = 0.5; % coherent proccessing interval (s)
        passive.cpi_overlap = 0.8; % overlap between CPIs (watch this - too large will cause slow exceution)
        passive.doppler_window = 'Blackman-Harris';
        passive.dopp_zero_padding = 1;
        passive.dynamic_range = +50;
        passive.range_bin_size = (1/passive.Fs)/passive.range_zero_padding * C/2;
        passive.max_range_m = passive.max_range*passive.range_bin_size;
        passive.min_range_m = passive.min_range*passive.range_bin_size;
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
        passive.range_axis = linspace(passive.min_range_m,passive.max_range_m,passive.no_range_bins);
        
        
        % create video of Passive range-Doppler slices        
        video_name = exp_dir + "passive_range-Doppler_" +'simple_running_ref_ideal' + ".avi";
        video_title = "Passive Pre-DSI";
        dynamic_range = 80;
        max_range = 100;
        max_doppler = 30;
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
        max_iterations = 1000; % maximum number of itterations DSI is CLEANed from CAF slice
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
        video_name = exp_dir + "CLEANed_range-Doppler_CLEANed_log_Exp_" +'simple_running' + ".avi";
        %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";
        video_title = "CLEANed Passive Radar Capture";
        dynamic_range = 80;
        max_range = 100;
        max_doppler = 30;
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



