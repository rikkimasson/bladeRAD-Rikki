clear all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 1;    % Expeiment Name
capture_duration = 1;        % capture duration
save_directory = "/home/piers/Documents/Captures/"; % rach experiment will save as a new folder in this directory


% FMCW Parameters 
FMCW_Fs = 40e6;          % Sample Rate of SDR per I & Q (in reality Fs is double this)
pulse_duration = 1e-3;   % Desired Pulse Duration 
FMCW_Bw = 40e6;          % LFM Bandwidth 
FMCW_Fc = 2400e6;   % Central RF 
Tx_gain = 20;       
Rx1_gain = 36;
Rx2_gain = 0;
Tx_SDR = 1;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
Rx_SDR = 2;   % SDR to use for RX
% Procesing Parameters
 max_range = 1000; %max range to LPF filter data to

% Passive Radar Parameters
passive_Fc = 5220e6;   % Central RF    
Ref_gain = 26; % 26dB seems good for C-Band Patch Antennas
Sur_gain = 40; % 40dB seems good for C-Band Patch Antennas
Pass_SDR = 3;   % SDR to use for Passive Radar - labelled on RFIC Cover and bladeRAD Facia Panel
passive_Bw = 20e6;
 
 
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


%% Create Sawtooth Chirp for bladeRF
chirp = saw_LFM_chirp(FMCW_Bw,pulse_duration,FMCW_Fs);
save_sc16q11('/tmp/chirp.sc16q11', chirp); %save chirp to binary file
clear chirp
    % spectrogram(chirp,128,100,128,Fs,'centered','yaxis') %plot spectrogram of chirp
    % f = linspace(-0.5 * Fs, 0.5 * Fs, length(chirp));
    % figure
    % plot(f,20*log10(abs(fftshift(fft(chirp)))/Samps_per_pulse));%plot FFT of chirp to show entire spectral content

%% Setup FMCW Radar
    % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''
 
    % Setup Tx SDR 
    [trig_flag_1,tx_command] = create_shell_command(Experiment_ID,...
                                   FMCW_number_cap_samps,... 
                                   FMCW_number_pulses,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   FMCW_RF_freq,...
                                   FMCW_Bw_M,...
                                   Tx_SDR,...
                                   'slave',...
                                   3,...
                                   'tx');
    tx_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    pause(5);    
    status = system(tx_command);


    % Setup Rx SDR 
    [trig_flag_2,rx_command] = create_shell_command(Experiment_ID,...
                                   FMCW_number_cap_samps,... 
                                   FMCW_number_pulses,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   FMCW_RF_freq,...
                                   FMCW_Bw_M,...
                                   Rx_SDR,...
                                   'slave',...
                                   1,...
                                   'rx'); 
    if trig_flag_1 && trig_flag_2
        "Trigger Conflict - FMCW Radar"
        return
    end
    rx_command = rx_command + "&"; % uncomment for non-blocking system command execution                                                              
    system(rx_command); % Blocking system command execution 
    pause(5);
    
%% Setup Passive Radar
    % by default clock out is always enabled 
    % 1 'set clock_sel external'; 
    % 2 'set clock_ref enable';

    % Setup Tx SDR 
   [trig_flag_3, passive_command] = create_shell_command(Experiment_ID,...
                                   passive_number_cap_samps,... 
                                   0,...
                                   0,...
                                   Ref_gain,...
                                   Sur_gain,...
                                   passive_RF_freq,...
                                   passive_Bw_M,...
                                   Pass_SDR,...
                                   'master',...
                                   1,...
                                   'pass');
    if trig_flag_1 && trig_flag_3 || trig_flag_2 && trig_flag_3
        "Trigger Conflict - Passive"
        return
    end
    %passive_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(passive_command);
    
    

%% Save Raw Data and create header to directory 
    % make save directory
        exp_dir = save_directory + Experiment_ID + '/';
        make_dir = 'mkdir ' + exp_dir;
        system(make_dir);
    % move FMCW receive file to save directory
        move_file = 'mv /tmp/fmcw_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
        rtn = system(move_file);
        if rtn == 0
            "FMCW Data Copyied to Save directory"
        else 
            "FMCW Copy Failed"
            return

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
    
% %% FMCW Processing and Print RTI
%     % load refsig for deramping
%         refsig = load_refsig(FMCW_Bw_M,FMCW_Fc,pulse_duration);    
%     % load Signal, Mix and Dermap Signal  
%         file_location = exp_dir + 'fmcw_' + Experiment_ID;
%         [max_range_actual,processed_signal] = deramp_and_decimate(file_location,max_range,refsig,capture_duration,FMCW_number_pulses,FMCW_Fs,slope);
%         save(exp_dir + 'deramped_signal','processed_signal')
%     % Plot RTI
%         Range_axis = linspace(0,max_range_actual,size(processed_signal,1));
%         Range_bin = 1:size(processed_signal,1);
%         time_axis = linspace(0,size(processed_signal,2)*pulse_duration,size(processed_signal,2));
%         RTI_plot= transpose(10*log10(abs(processed_signal./max(processed_signal(:)))));
%         figure
%         imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
        
 %% Passsive Processing
    %load signal and split ref and sur
        file_location = exp_dir + 'passive_' + Experiment_ID;
        [ref_channel, sur_channel]  = load_passive_data(file_location);
    % Plot time domain signals
         figure
         plot(real(ref_channel))
         title("Ref channel time series");
         figure
         plot(real(sur_channel))
         title("Sur channel time series");       