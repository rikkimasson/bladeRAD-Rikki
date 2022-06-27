clear all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 22;    % Expeiment Name
capture_duration = 15;        % capture duration
% save_directory = "/media/piers/data_drive/BladeRF_Experiments/Hybrid Radar/"; % each experiment will save as a new folder in this directory
save_directory = "/home/piers/Documents/Captures/7_Oct/hybrid/"; % rach experiment will save as a new folder in this directory

Experiment_ID = 2;    % Expeiment Name
capture_duration = 15;        % capture duration
% save_directory = "/media/piers/data_drive/BladeRF_Experiments/Hybrid Radar/"; % each experiment will save as a new folder in this directory
% save_directory = "/home/piers/Documents/Captures/6_Oct/hybrid/"; % rach experiment will save as a new folder in this directory
save_directory = "/media/piers/data_drive/BladeRF_Experiments/6_Oct/hybrid/"

exp_dir = save_directory + Experiment_ID + '/';


% FMCW Parameters 
FMCW_Fs = 20e6;          % Sample Rate of SDR per I & Q (in reality Fs is double this)
pulse_duration = 1e-3;   % Desired Pulse Duration 
FMCW_Bw = 20e6;          % LFM Bandwidth 
FMCW_Fc = 5.8e9;   % Central RF 
Tx_gain = 66;       % [-23.75, 66] (S-Band = 23.5 dBm) (C-Band = 15.8 dBm)
Rx1_gain = 50;      % [-16, 60]
Rx2_gain = 0;       % [-16, 60]
Rx_1_lna = True;
Rx_2_lna = True;
Tx_SDR = 1;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
Rx_SDR = 2;   % SDR to use for RX
% Procesing Parameters
FMCW_max_range = 1000; %max range to LPF filter data to

% Passive Radar Parameters
passive_Fc = 2422e6;   % Central RF    
Ref_gain = 5; % 26dB seems good for C-Band Patch Antennas [-16, 60]
Sur_gain = 30; % 40dB seems good for C-Band Patch Antennas [-16, 60]
Pass_SDR = 3;   % SDR to use for Passive Radar - labelled on RFIC Cover and bladeRAD Facia Panel
passive_Bw = 20e6;
passive_max_range = 100; %max range to cross-correlate to
 
% Parameters not configurable by user 
    C = physconst('LightSpeed');
    % FMCW parameters
    PRF = 1/pulse_duration;
    slope = FMCW_Bw/pulse_duration;
    F_Max = FMCW_Bw/2;
    R_Max = beat2range(F_Max,slope);
    FMCW_sample_duration = 1/FMCW_Fs;
    samples_per_pulse = pulse_duration/FMCW_sample_duration;
    FMCW_number_pulses = capture_duration / pulse_duration;
    FMCW_number_cap_samps = capture_duration/FMCW_sample_duration;
    FMCW_RF_freq = FMCW_Fc/1e6;   % RF in MHz 
    FMCW_Bw_M = FMCW_Bw/1e6;      % BW in MHz
    % passive parameters
    passive_Fs = passive_Bw;
    passive_sample_duration = 1/passive_Fs;
    passive_number_cap_samps = (2*capture_duration)/passive_sample_duration;
    passive_RF_freq = passive_Fc/1e6;   % RF in MHz 
    passive_Bw_M = passive_Bw/1e6;      % BW in MHz


% %% Create Sawtooth Chirp for bladeRF
% chirp = saw_LFM_chirp(FMCW_Bw,pulse_duration,FMCW_Fs);
% save_sc16q11('/tmp/chirp.sc16q11', chirp); %save chirp to binary file
% clear chirp
% 
% %% Setup FMCW Radar
%     % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''
%     
%     %initialise_clk_dist.sh
%  
%     % Setup Tx SDR 
%     [trig_flag_1,tx_command] = create_shell_command(Experiment_ID,...
%                                    FMCW_number_cap_samps,... 
%                                    FMCW_number_pulses,...
%                                    0,...
%                                    Tx_gain,...
%                                    Rx1_gain,...
%                                    Rx2_gain,...
%                                    FMCW_RF_freq,...
%                                    FMCW_Bw_M,...
%                                    Tx_SDR,...
%                                    'slave',...
%                                    3,...
%                                    'tx');
%     tx_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
%     status = system(tx_command);
%     pause(5);
% 
% 
%     % Setup Rx SDR 
%     [trig_flag_2,rx_command] = create_shell_command(Experiment_ID,...
%                                    FMCW_number_cap_samps,... 
%                                    FMCW_number_pulses,...
%                                    0,...
%                                    Tx_gain,...
%                                    Rx1_gain,...
%                                    Rx2_gain,...
%                                    FMCW_RF_freq,...
%                                    FMCW_Bw_M,...
%                                    Rx_SDR,...
%                                    'slave',...
%                                    1,...
%                                    'rx'); 
%     if trig_flag_1 && trig_flag_2
%         "Trigger Conflict - FMCW Radar"
%         return
%     end
%     rx_command = rx_command + "&"; % uncomment for non-blocking system command execution                                                              
%     system(rx_command); % Blocking system command execution 
%     pause(5);
%     
% %% Setup Passive Radar
%     % by default clock out is always enabled 
%     % 1 'set clock_sel external'; 
%     % 2 'set clock_ref enable';
% 
%     % Setup Passive SDR 
%    [trig_flag_3, passive_command] = create_shell_command(Experiment_ID,...
%                                    passive_number_cap_samps,... 
%                                    0,...
%                                    0,...
%                                    0,...
%                                    Ref_gain,...
%                                    Sur_gain,...
%                                    passive_RF_freq,...
%                                    passive_Bw_M,...
%                                    Pass_SDR,...
%                                    'master',...
%                                    1,...
%                                    'pass');
%                                
%     if trig_flag_1 && trig_flag_3 || trig_flag_2 && trig_flag_3
%         "Trigger Conflict - Passive"
%         return
%     end
%     %passive_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
%     status = system(passive_command);
%     
%     
% 
% %% Save Raw Data and create header to directory 
%     % make save directory
%         make_dir = 'mkdir ' + exp_dir;
%         system(make_dir);
%     % move FMCW receive file to save directory
%         move_file = 'mv /tmp/active_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
%         rtn = system(move_file);
%         if rtn == 0
%             "FMCW Data Copyied to Save directory"
%         else 
%             "FMCW Copy Failed"
%             return
% 
%         end
%     % move passive file to save directory
%     move_file = 'mv /tmp/passive_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
%     rtn = system(move_file);
%     if rtn == 0
%         "Passive Data Copyied to Save directory"
%     else 
%         "Passive Copy Failed"
%         return
%     end
%     save(exp_dir + 'Experimental Configuration') 
%     
%% FMCW Processing and Print RTI
    % load refsig for deramping
        refsig = load_refsig(FMCW_Bw_M,FMCW_Fc,pulse_duration);    
    % load Signal, Mix and Dermap Signal  
        file_location = exp_dir + 'active_' + Experiment_ID;
        zero_padding = 2; % 1 = none; 2 = 100%
        [max_range_actual,processed_signal] = deramp_and_decimate(file_location,FMCW_max_range,refsig,capture_duration,FMCW_number_pulses,FMCW_Fs,slope,zero_padding);
        save(exp_dir + 'deramped_signal','processed_signal')
    % Plot RTI
        Range_axis = linspace(0,max_range_actual,size(processed_signal,1));
        Range_bin = 1:size(processed_signal,1);
        time_axis = linspace(0,size(processed_signal,2)*pulse_duration,size(processed_signal,2));
        RTI_plot= transpose(10*log10(abs(processed_signal./max(processed_signal(:)))));
        figure
        fig = imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
            ylabel('Time (Sec)')
            xlabel('Range Bin')
            title("FMCW RTI - " + Experiment_ID)
            xlim([0 100])
            fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)    


      % Spectrogram 
        r_bin = 1;
        l_fft = 2056;
        pad_factor = 4;
        overlap_factor = 0.99;
        [spect,f] = spectrogram(processed_signal(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/passive_Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-30 0]);   
            ylim([-100 100])
            colorbar
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "FMCW Spectrogram - R Bin: " + r_bin + " - " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)


 %% Passive Processing
    % load signal and split ref and sur
        file_location = exp_dir + 'passive_' + Experiment_ID;
        [ref_channel, sur_channel]  = load_passive_data(file_location);
    % Plot time domain signals
         figure
         fig = subplot(2,1,1);
            plot(real(ref_channel(1:400000)));
            title("Ref channel time series");
         hold on
         subplot(2,1,2)
             plot(real(sur_channel(1:400000)))
             title("Sur channel time series");    
             fig_name = exp_dir + "Time Domain Signals_" + Experiment_ID + ".jpg";
             saveas(fig,fig_name,'jpeg')
   % Batch process data and cross correlate  
         seg_s = 5000; % number of segments per second - analagos to PRF.
         seg_percent = 10;  % percentage of segment used for cross coreclation of 
                            % survallance and reference. Will affect SNR dramatically.
         cc_matrix = passive_batch_process(ref_channel,sur_channel,seg_s,seg_percent,passive_Fs,passive_max_range,exp_dir);
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
        

<<<<<<< HEAD
        end
    % move passive file to save directory
    move_file = 'mv /tmp/passive_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Passive Data Copyied to Save directory"
    else 
        "Passive Copy Failed"
        return
    end
    save(exp_dir + 'Hybrid Experimental Configuration') 
    
%% FMCW Processing and Print RTI
    % load refsig for deramping
        refsig = load_refsig(FMCW_Bw_M,FMCW_Fc,pulse_duration);    
    % load Signal, Mix and Dermap Signal  
        file_location = exp_dir + 'active_' + Experiment_ID;
        zero_padding = 2; % 1 = none; 2 = 100%
        [max_range_actual,processed_signal] = deramp_and_decimate(file_location,FMCW_max_range,refsig,capture_duration,FMCW_number_pulses,FMCW_Fs,slope,zero_padding);
        save(exp_dir + 'deramped_signal','processed_signal')
    % Plot RTI
        Range_axis = linspace(0,max_range_actual,size(processed_signal,1));
        Range_bin = 1:size(processed_signal,1);
        time_axis = linspace(0,size(processed_signal,2)*pulse_duration,size(processed_signal,2));
        RTI_plot= transpose(10*log10(abs(processed_signal./max(processed_signal(:)))));
        figure
        fig = imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
            ylabel('Time (Sec)')
            xlabel('Range Bin')
            title("FMCW RTI - " + Experiment_ID)
            xlim([0 100])
            fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)    


      % Spectrogram 
        start_bin = 3;
        finish_bin = 9;
        spect_wind = sum(processed_signal(start_bin:finish_bin,:));
        l_fft = 512;
        pad_factor = 4;
        overlap_factor = 0.99;
        [spect,f] = spectrogram(spect_wind,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/passive_Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-30 0]);   
            ylim([-100 100])
            colorbar
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "FMCW Spectrogram - " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)


 %% Passive Processing
    % load signal and split ref and sur
        file_location = exp_dir + 'passive_' + Experiment_ID;
        [ref_channel, sur_channel]  = load_passive_data(file_location);
    % Plot time domain signals
         figure
         fig = subplot(2,1,1);
            plot(real(ref_channel(1:400000)));
            title("Ref channel time series");
         hold on
         subplot(2,1,2)
             plot(real(sur_channel(1:400000)))
             title("Sur channel time series");    
             fig_name = exp_dir + "Time Domain Signals_" + Experiment_ID + ".jpg";
             saveas(fig,fig_name,'jpeg')
   % Batch process data and cross correlate  
         seg_s = 5000; % number of segments per second - analagos to PRF.
         seg_percent = 10;  % percentage of segment used for cross coreclation of 
                            % survallance and reference. Will affect SNR dramatically.
         cc_matrix = passive_batch_process(ref_channel,sur_channel,seg_s,seg_percent,passive_Fs,passive_max_range,exp_dir);
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
        start_bin = 1;
        finish_bin = 6;
        spect_wind = sum(cc_matrix(start_bin:finish_bin,:));
        l_fft = 512;
        pad_factor = 4;
        overlap_factor = 0.99;
        [spect,f] = spectrogram(spect_wind,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,seg_s,'centered','yaxis');
=======
     % Spectrogram 
        r_bin = 1;
        l_fft = 2056;
        pad_factor = 4;
        overlap_factor = 0.99;
        [spect,f] = spectrogram(cc_matrix(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,seg_s,'centered','yaxis');
>>>>>>> 1c5ddb33d1ae92180b94ba0050155e51f243fe60
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/passive_Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-30 0]);   
<<<<<<< HEAD
            ylim([-100 100])
=======
            ylim([-600 600])
>>>>>>> 1c5ddb33d1ae92180b94ba0050155e51f243fe60
            colorbar
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
<<<<<<< HEAD
            fig_title = "Passive Spectrogram -  "  + Experiment_ID;
=======
            fig_title = "Passive Spectrogram - R Bin: " + r_bin + " - " + Experiment_ID;
>>>>>>> 1c5ddb33d1ae92180b94ba0050155e51f243fe60
            title(fig_title);
            fig_name = exp_dir + "Passive Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)



         
