
% 1. Respond in the command window when prompted in order to point MATLAB
% to the correct directory
% prompt = "Which computer is running the script? 1=Rikki's Mac, 2=Jacob's Laptop, 3=Rikki UCL \n"

clc
clear
close all
tic
pc = 3;%input(prompt);
addpath("C:/Users/rw/Documents/RADCOM/misc_helper_functions");
if pc == 1
    local_save_directory = "/Users/rikkimasson/Documents/MATLAB/experiment_data/";
    repo_directory = "/Users/rikkimasson/Documents/GitHub";
    elseif pc == 2
        local_save_directory = "C:\Users\jster\OneDrive - University College London\Processing Project Jacob + Rikki\Data Files";
        repo_directory = "D:\GitHub Repos";
    elseif pc == 3
        local_save_directory = "D:\Passive_Data\Data Files";
        repo_directory = "C:\Users\rw\Documents\";
     elseif pc == 4
        local_save_directory = "C:\Users\uceerm6\OneDrive - University College London\Documents\Radar_Data\beasley_bladeRAD\";
        repo_directory = "C:\Users\uceerm6\OneDrive - University College London\Documents\GitHub";
    else error('Error. \nInput number must correspond to an option') 
end

% 2. Processing Flags - choose whether to
process_active_flag = false;
process_passive_flag = true;
ground_truth_flag = false;

% 3. Select the experiment you wish to process
experiment_number = 2;

if ismac
        bracket = "/";
    elseif isunix
        disp('isunix repo directory not defined')
    elseif ispc
        bracket = "\";
    else
        disp('Platform not supported')
end

% adds paths to generic functions
if ismac
    addpath(repo_directory + '/bladeRAD-Rikki/generic_scripts/matlab',...
    repo_directory + '/bladeRAD-Rikki/generic_scripts/matlab/CFAR/',...
    repo_directory + '/bladeRAD-Rikki/generic_scripts',...
    repo_directory + '/bladeRAD-Rikki/generic_scripts/ref_signals/')
    elseif isunix
            disp('isunix repo directory not defined')
    elseif ispc
        addpath(repo_directory + '\bladeRAD-Rikki\generic_scripts\matlab',...
        repo_directory + '\bladeRAD-Rikki\generic_scripts\matlab\CFAR\',...
        repo_directory + '\bladeRAD-Rikki\generic_scripts',...
        repo_directory + '\bladeRAD-Rikki\generic_scripts\ref_signals\')
    else
        disp('Platform not supported')
end

for i=experiment_number
%% Load .mat file containing experiment setup parameters

    mat_file_name = local_save_directory + i + bracket + "Experimental Configuration.mat";
    load(mat_file_name);
    exp_dir = local_save_directory + i + bracket;
%% Load Ground Truth Data

    if ground_truth_flag == true
        mat_file_name = local_save_directory  + i + bracket +"geode1.mat";
        load(mat_file_name);
        
        gt.time_offset = 280-14; % offset between start time of radar and gps recording
        gt.range_offset = 8; % Range fudge factor between radar and gps
        gt.range = geode1.target_range(gt.time_offset:gt.time_offset+capture_duration) + gt.range_offset;
        gt.sog = geode1.SoG(gt.time_offset :gt.time_offset +capture_duration);
        
        % Interpolate to number of active radar pulses
        old_sample_points = 0:capture_duration;
        active_pulses_sample_points = 1:active.number_pulses;
        gt.pulses_interp_range = interp1(old_sample_points,gt.range,active_pulses_sample_points/1000,'linear');
    end
%% FMCW Processing and Print RTI

    if process_active_flag == true
        % load refsig for deramping
        refsig = load_refsig(active.Bw_M,active.Fs,active.Fc,active.pulse_duration);
        % load Signal, Mix and Dermap Signal
        file_location = exp_dir + 'active_' + Experiment_ID;
        lp_filter = true;
        [max_range_actual,deramped_signal,active.decimation_factor_actual] = deramp_and_decimate(file_location,active.max_range,refsig,capture_duration,active.number_pulses,active.Fs,active.slope,lp_filter);
        save(exp_dir + 'deramped_signal','deramped_signal')
        
        % Window and FFT Signal
        % window signal
        w = window('blackman',size(deramped_signal,1));
        windowed_signal = deramped_signal.*w;
        % fft signal
        zero_padding = 1; % 1 = none; 2 = 100%
        processed_signal = fft(windowed_signal,size(windowed_signal,1)*zero_padding);
        beat_frequncies = processed_signal(1:(size(processed_signal,1)/2),:); % keep +ve beat frequencies
        
        % MTI Filtering
        % Single Delay Line Filter
        MTI_Data = zeros(size(beat_frequncies));
        active.range_bins = size(MTI_Data,1);
        for i=2:active.number_pulses
            MTI_Data(:,i) = beat_frequncies(:,i)-beat_frequncies(:,i-1);
        end
        % IIR Filter
        [b, a] = butter(10, 0.021, 'high');
        for i=1:active.range_bins
            MTI_Data(i,:) = filtfilt(b,a,beat_frequncies(i,:));
        end
        
        % Derive range and time axis
        active.n_range_bins = size(beat_frequncies,1);
        active.range_bins = 1:active.n_range_bins;
        active.fftfrequncies =fftfreq(size(processed_signal,1),1/(active.Fs/active.decimation_factor_actual)); % possible beat frequencies
        active.slope = active.Bw/active.pulse_duration;
        ranges = (active.fftfrequncies*C)/(2*active.slope); % calculate true range bin size
        active.range_axis = ranges(1:(size(ranges,2)/2)); % truncate to only +ve beat frequencies
        active.range_bin_size = ranges(2)
        active.time_axis = linspace(0,size(processed_signal,2)*active.pulse_duration,size(processed_signal,2));
        
        % Plot RTI
        RTI_plot= transpose(20*log10(abs(beat_frequncies./max(beat_frequncies(:)))));
        figure
        fig = imagesc(active.range_axis,active.time_axis,RTI_plot,[-50,0]); hold on;
        if ground_truth_flag == true
            plot(gt.pulses_interp_range,active.time_axis,'-r', 'LineWidth', 2);
        end
        ylabel('Time (s)');xlabel('Range (m)');legend('Target Ground Truth') % Ground Truth data not provided
        title("FMCW RTI - " + Experiment_ID)
        xlim([0 200])
        fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')
        saveas(fig,fig_name)
        
        % Plot MTI RTI
        RTI_plot= transpose(20*log10(abs(MTI_Data./max(MTI_Data(:)))));
        figure
        fig = imagesc(active.range_axis,active.time_axis,RTI_plot,[-50,0]); hold on;
        if ground_truth_flag == true
            plot(gt.pulses_interp_range,active.time_axis,'-r', 'LineWidth', 2);
        end
        ylabel('Time (s)');xlabel('Range (m)');legend('Target Ground Truth') % Ground Truth data not provided
        title("FMCW RTI - " + Experiment_ID)
        xlim([0 200])
        fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')
        saveas(fig,fig_name)
        
        % Plot series of pulses
        fig = figure
        ranges_2_plot = floor(linspace(1,active.number_pulses,5));
        for i = ranges_2_plot
            plot(active.range_axis,RTI_plot(i,:));
            hold on
        end
        title("Single Pulse - " + Experiment_ID);
        xlim([0 200])
        grid on; grid minor;
        ylabel('Relative Power (dB)')
        xlabel('Range (m)')
        fig_name = exp_dir + "Single_Pulse" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')
        
        
        % Spectrogram
        % Parameters
        r_start = 5;
        r_stop = 12;
        l_fft = 256;
        pad_factor = 1;
        overlap_factor = 0.90;
        
        % Plot Spectrogram pre-MTI filtering
        integrated_data = sum(beat_frequncies(r_start:r_stop,:));
        [spect,active.doppler_axis] = spectrogram(integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,active.PRF,'centered','yaxis');
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(active.time_axis,-active.doppler_axis,spect,[-50 0]);
        ylim([-500 500])
        c = colorbar;
        c.Label.String='Norm Power (dB)'
        xlabel('Time (Sec)')
        % ylabel('Radial Velocity (mph)')
        ylabel('Doppler Frequency (Hz)')
        fig_title = "FMCW Spectrogram -" + Experiment_ID;
        title(fig_title);
        fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')
        fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID;
        saveas(fig,fig_name)
        
        % Plot Spectrogram post-MTI filtering
        MTI_integrated_data = sum(MTI_Data(r_start:r_stop,:));
        [spect,f] = spectrogram(MTI_integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,active.PRF,'centered','yaxis');
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(active.time_axis,-f,spect,[-50 0]);
        ylim([-500 500])
        c = colorbar
        c.Label.String='Norm Power (dB)'
        xlabel('Time (Sec)')
        % ylabel('Radial Velocity (mph)')
        ylabel('Doppler Frequency (Hz)')
        fig_title = "MTI FMCW Spectrogram - " + Experiment_ID;
        title(fig_title);
        fig_name = exp_dir + "MTI FMCW Spectrogram_" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')
        fig_name = exp_dir + "MTI FMCW Spectrogram_" + Experiment_ID;
        saveas(fig,fig_name)
%% Process Active data into Range-Doppler Slices

        active.cpi = 0.5; % cohernet proccessing interval (s)
        active.cpi_overlap = 0.8; % overlap between CPIs (watch this - too large will cause slow exceution)
        active.zero_padding = 1; % 1 = none; 2 = 100%
        active.doppler_window = 'Hann';
        active.dynamic_range = 50
        active.velocity_conv = C*(((1/C)/(active.Fc/C)));
        active.max_range_actual = max_range_actual; %1000;
        
        % Create range-Doppler surfaces
        [active.number_cpi,....
            active.pulses_per_cpi,...
            active.range_doppler_slices] = rangeDopplerSlice(MTI_Data,active.cpi,active.PRF,...
            active.cpi_overlap,...
            active.zero_padding,...
            active.doppler_window);
        
        % Create Active Data Range and Doppler Axis
        active.doppler_bins = 1:size(active.range_doppler_slices{1},1);
        active.doppler_axis = linspace(-active.PRF/2,active.PRF/2,size(active.doppler_bins,2));
        active.doppler_velocity_axis = active.doppler_axis * active.velocity_conv;
        active.range_bins_axis = 1:active.range_bins;
        
        % Create video of range-Doppler slices
        video_name = exp_dir + "active_range-Doppler_" + Experiment_ID + ".avi";
        % video_name = "range-Doppler_log_Exp_" + Experiment_ID + ".avi";
        video_title = "Active Radar Capture";
        dynamic_range = +100;
        max_range = 200;
        max_doppler = 60;
        frame_rate = 1/(capture_duration/active.number_cpi);
        createVideo(active.range_doppler_slices,frame_rate,...
            active.range_axis-10,max_range,...
            -active.doppler_velocity_axis,max_doppler,...
            dynamic_range,video_name,video_title);
        
    end
%% Passive Processing

    if process_passive_flag == true
        passive.max_range = 50; % maximum number of range bins - xcorr shifts.
        passive.range_zero_padding = 1; % 1 = none, 2 = 100%
        passive.td_corr = false; % true = time domain xcorr; false = freq domain xcorr
        passive.seg_s = 1000; % number of segments per second - analagos to PRF.
        passive.seg_percent = 100; % % of segment used for xcorr.
        
        % load signal and split ref and sur
        file_location = exp_dir + 'passive_' + Experiment_ID;
        [ref_channel, sur_channel]  = load_passive_data(file_location);
        % Plot time domain signals
        % figure
        % fig = subplot(2,1,1);
        % plot(real(ref_channel(1:4000000)));
        % title("Ref channel time series");
        % hold on
        % subplot(2,1,2)
        % plot(real(sur_channel(1:4000000)))
        % title("Sur channel time series");
        % fig_name = exp_dir + "Time Domain Signals_" + Experiment_ID + ".jpg";
        % saveas(fig,fig_name,'jpeg')
        
        % Add Noise to Passive Radar
        if (0)
        noisy_sur_channel= awgn(sur_channel,-2,'measured');
        else
            noisy_sur_channel=sur_channel;
        end
        % 
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
        


        % [ref_matrix_fpga ,self_ambg_matrix_fpga, cc_matrix_fpga] = passive_batch_process_fpga(ref_channel,noisy_sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.max_range,exp_dir,passive.range_zero_padding,passive.td_corr);
       

        % self_ambg_matrix=self_ambg_matrix_fpga;
        % cc_matrix=cc_matrix_fpga;
        % ref_matrix= ref_matrix_fpga;

        [ref_matrix ,self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,noisy_sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.max_range,exp_dir,passive.range_zero_padding,passive.td_corr);
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
        video_name = exp_dir + "passive_range-Doppler_" + Experiment_ID + ".avi";
        video_title = "Passive Pre-DSI";
        dynamic_range = +inf;
        max_range = 200;
        max_doppler = 100;
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
        video_name = exp_dir + "CLEANed_range-Doppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";
        %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";
        video_title = "CLEANed Passive Radar Capture";
        dynamic_range = 50;
        max_range = 100;
        max_doppler = 20;
        frame_rate = 1/(capture_duration/passive.number_cpi);
        createVideo(passive.CLEANed_range_doppler_slices,frame_rate,...
            passive.range_axis,max_range,...
            passive.doppler_velocity_axis,max_doppler,...
            dynamic_range,video_name,video_title);        
    end 
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



