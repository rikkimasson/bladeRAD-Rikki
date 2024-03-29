clear all
close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 13;       % Expeiment Name
capture_duration = 30;    % capture duration
% save_directory = "/media/sdrlaptop1/T7/22_06_21_N0/"; % each experiment will save as a new folder in this directory
save_directory = "~/Documents/bladeRAD_Captures/06_07_2022_farm/hybrid_radar/"; % each experiment will save as a new folder in this directory

passive.max_range = 100; %max range to cross-correlate to

% Radar Parameters 
passive.Fc = 2437e6;   % Central RF
passive.Bw = 20e6;               % Sample Rate of SDR per I & Q (in reality Fs is double this)
passive.Ref_gain = 3;
passive.Sur_gain = 3;
passive.SDR = 3;   % SDR to use for Passive Radar - labelled on RFIC Cover and bladeRAD Facia Panel

% Parameters not configurable by user 
    C = physconst('LightSpeed');
    passive.Fs = passive.Bw;
    passive.sample_duration = 1/passive.Fs;
    passive.number_cap_samps = 2*(capture_duration/passive.sample_duration)
    passive.RF_freq = passive.Fc/1e6;   % RF in MHz 
    passive.Bw_M = passive.Bw/1e6;      % BW in MHz

    passive_file_size_MBytes = (passive.number_cap_samps * 16)*2/(8*1e6) 

    
%% Setup Radar
    % 1 'set clock_sel external'; 2 'set clock_out enable; 3 'set clock_ref enable'

    % Setup Passive SDR 
   [trig_flag, passive_command] = create_shell_command(Experiment_ID,...
                                   passive.number_cap_samps,... 
                                   0,...
                                   0,...
                                   0,...
                                   passive.Ref_gain,...
                                   passive.Sur_gain,...
                                   passive.RF_freq,...
                                   passive.Bw_M,...
                                   passive.SDR,...
                                   'master',...
                                   3,...
                                   'pass');
                               
    %passive_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(passive_command);

%% Save Raw Data and create header file for directory 
    exp_dir = save_directory + Experiment_ID + '/';
    make_dir = 'mkdir ' + exp_dir;
    system(make_dir); % Blocking system command execution
    move_file = 'mv /tmp/passive_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Rx Data Copyied to Save directory"
    else 
        "Rx Copy Failed"
        return
    end
    save(exp_dir + 'Passive Experimental Configuration') 


 %% Passive Processing
    % load signal and split ref and sur
        file_location = exp_dir + 'passive_' + Experiment_ID;
        [ref_channel, sur_channel]  = load_passive_data(file_location);
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
         passive.seg_s = 5000; % number of segments per second - analagos to PRF.
         passive.seg_percent = 100;  % percentage of segment used for cross coreclation of 
                            % survallance and reference. Will affect SNR dramatically.
         [ref_matrix ,self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.max_range,exp_dir);
         save(exp_dir + 'passive_matrix','cc_matrix')
    % RTI Plot
        RTI_plot= transpose(10*log10(abs(cc_matrix./max(cc_matrix(:)))));
        Range_bin = linspace(0,passive.max_range,size(cc_matrix,1));
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
        f_axis = linspace(-passive.seg_s/2,passive.seg_s/2,size(cc_matrix,2));
        t_cc_matrix = transpose(cc_matrix);
        CAF = fftshift(fft(t_cc_matrix,size(t_cc_matrix,1),1),1);
        figure
        imagesc(Range_bin,f_axis,10*log10(abs(CAF./max(CAF(:)))),[-50 1]); 
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
        r_start = 1;
        r_stop = 100;
        l_fft = 1024;
        pad_factor = 4;
        overlap_factor = 0.99;
        integrated_data = sum(t_cc_matrix(r_start:r_stop,:));
        [spect,f] = spectrogram(cc_matrix(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,passive.seg_s,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/passive.Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-50 0]);   
            ylim([-50 50])
            colorbar
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "Passive Spectrogram - R Bin: " + r_bin + " - " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "Passive Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)


%% Proccess Passive data into Range-Doppler Slices
           passive.PRF = passive.seg_s; %seg_s
           passive.cpi = 0.5; % cohernet proccessing interval (s)
           passive.cpi_overlap = 0.5; % overlap between CPIs (watch this - too large will cause slow exceution)
           passive.doppler_window = 'hann';
           passive.zero_padding = 4;
           passive.dynamic_range = +inf;
           passive.max_range = passive.max_range;       
     
                [passive.number_cpi,... 
                 passive.pulses_per_cpi,...
                 passive.range_doppler_slices] = rangeDopplerSlice(cc_matrix,passive.cpi,passive.PRF,...
                                                                   passive.cpi_overlap,...
                                                                   passive.zero_padding,...
                                                                   passive.doppler_window);  
                
                 [~,~,passive.self_ambg_slices] = rangeDopplerSlice(self_ambg_matrix,passive.cpi,passive.PRF,...
                                                                   passive.cpi_overlap,...
                                                                   passive.zero_padding,...
                                                                   passive.doppler_window);  
            
            passive.cpi_stride = round(passive.pulses_per_cpi*(1-passive.cpi_overlap)); % number of pulses to stride each for next CPI
            passive.velocity_conv = C*(((1/C)/(passive.Fc/C)));
            passive.range_bins = size(cc_matrix,1);
            passive.doppler_bins = passive.pulses_per_cpi*passive.zero_padding+1;
            passive.doppler_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.doppler_bins);
            passive.doppler_velocity_axis = passive.doppler_axis*passive.velocity_conv;


            % create video of CLEANed range-Doppler slices
             video_name = exp_dir + "Range_Doppler_Slices" + Experiment_ID + ".avi";
             %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";       
             video_title = "Passive Pre-DSI";
             dynamic_range = 50;
             max_range = 50;
             max_doppler = 50;
             frame_rate = 1/(capture_duration/passive.number_cpi);    
             createVideo(passive.range_doppler_slices,frame_rate,...
                         passive.range_axis,max_range,...
                         passive.doppler_axis,max_doppler,...
                         dynamic_range,video_name,video_title);
        
                 
%% Direct Signal Interference Cancellation
   % set DSI cancellation parameters
      p = 0.999;        % subtraction parameter - P must be a positive 
                        % integer less than one to avoid unwanted discontinuities
                        % arising from zero values in the range–Doppler surface.
     threshold = 0.005; % cutoff threshold parameter
     max_iterations = 20; % maximum number of itterations DSI is CLEANed from CAF slice
     number_rbins = size(passive.range_doppler_slices{1},2);
  
   % perform CLEAN based DSI Cancellation   
     passive.CLEANed_range_doppler_slices = CLEAN_DSI(passive.range_doppler_slices,...
                                               passive.self_ambg_slices,...
                                               ref_matrix,...
                                               passive.number_cpi,...
                                               passive.pulses_per_cpi,...
                                               passive.cpi_stride,...
                                               passive.zero_padding,...
                                               number_rbins,...
                                               max_iterations,...
                                               threshold,p,...
                                               passive.range_axis,passive.doppler_axis);


  % create video of CLEANed range-Doppler slices
     video_name = exp_dir + "CLEANed_range-Doppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";
     %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";       
     video_title = "CLEANing Passive Radar Capture";
     dynamic_range = 50;
     max_range = 50;
     max_doppler = 500;
     frame_rate = 1/(capture_duration/passive.number_cpi);    
     createVideo(passive.CLEANed_range_doppler_slices,frame_rate,...
                 passive.range_axis,max_range,...
                 passive.doppler_axis,max_doppler,...
                 dynamic_range,video_name,video_title);
