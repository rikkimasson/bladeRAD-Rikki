% clear all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 3;    % Expeiment Name
capture_duration = 20;        % capture duration
save_directory = "/media/piers/data_drive/BladeRF_Experiments/Hybrid Radar/"; % each experiment will save as a new folder in this directory
% save_directory = "/home/piers/Documents/Captures/"; % rach experiment will save as a new folder in this directory

% Radar Parameters
Fc = 2440e6;
Fs = 40e6;          %Sample Rate of SDR per I & Q (in reality Fs is double this)
pulse_duration = 100e-6;   %Desired Pulse Duration 
Bw = 40e6;          %LFM Bandwidth 
PRF = 1000;
Tx_gain = 20;       
Rx1_gain = 36;
Rx2_gain = 0;
Tx_SDR = 1;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
Rx_SDR = 2;   % SDR to use for RX

% Parameters not configurable by user 
    C = physconst('LightSpeed');
    PRI = 1/PRF;
    if PRI > pulse-duration
        "pulse duration longer than PRI"
        return
    end
    number_pulses = capture_duration/PRI;
    sample_duration = 1/Fs;
    samples_per_pulse = pulse_duration/sample_duration;
    samples_per_PRI = PRI/sample_duration;
    Num_cap_samps = capture_duration/sample_duration;
    Tx_Delay = PRI - pulse_duration % in seconds
    Tx_Delay_us = Tx_Delay * 10e5 % in micro seconds
    R_Max = PRI*C/2;
    Fc_M = Fc/1e6;      % RF in MHz 
    Bw_M = Bw/1e6;      % BW in MHz

%% Create Sawtooth Chirp for bladeRF
 chirp = saw_LFM_chirp(Bw,pulse_duration,Fs);
 save_sc16q11('/tmp/chirp.sc16q11', chirp); %save chirp to binary file
 clear chirp


%% Setup Radar
    % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''
 
    % Setup Tx SDR 
    [trig_flag_1,tx_command] = create_shell_command(Experiment_ID,...
                                   0,... 
                                   number_pulses,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   Fc_M,...
                                   Bw_M,...
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
                                   0,...
                                   0,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   Fc_M,...
                                   Bw_M,...
                                   Rx_SDR,...
                                   'slave',...
                                   1,...
                                   'rx'); 
    if trig_flag_1 && trig_flag_2
        "Trigger Conflict"
        return
    end
    rx_command = rx_command + "&"; % uncomment for non-blocking system command execution                                                              
    system(rx_command); % Blocking system command execution 
    pause(5);
 
%% Save Raw Data and create  header to directory 
    exp_dir = save_directory + Experiment_ID + '/';
    make_dir = 'mkdir ' + exp_dir;
    system(make_dir); % Blocking system command execution
    move_file = 'mv /tmp/fmcw_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Rx Data Copyied to Save directory"
    else 
        "Rx Copy Failed"
        return
        
    end
    save(exp_dir + 'FMCW Experimental Configuration') 

  %% Load Reference Deramp Signal
    refsig = load_refsig(Bw_M,Fc,pulse_duration);
   
%% Load Signals
    file_location = exp_dir + 'active_' + Experiment_ID;
    raw_data = load_sc16q11(file_location);

    
%% Reshape array into matrix of pulses
pulse_matrix = reshape(raw_data,[length(raw_data)/number_pulses,number_pulses]); %reshape array to individual pulses

spectrogram(raw_data(1:4000),128,100,100,Fs,'centered','yaxis')

%% Range Limmited Signal Processing 

zero_padding = 1;

% Match Filter 
matched_filter = transpose(refsig);
% refsig = refsig(1:samples_per_pulse);


matched_filter_fft = conj(fft(matched_filter,(samples_per_PRI*zero_padding))); 
radar_matrix = zeros(size(pulse_matrix,1),size(matched_filter_fft,1));

for i = 1:number_pulses
    appo = fft(pulse_matrix(1:end,i));
    radar_matrix(:,i) = ifft(appo.*matched_filter_fft); 
end

%Range_bin = 1:size(Data_matched,1) * samp_dur * C;
    Range_bin = 1:size(radar_matrix,1);
    Range = 1:R_Max;
    time_axis = 1:capture_duration; %linspace(0,PRI,size(N0_Data_matched,2));
    RTI_plot= transpose(10*log10(abs(radar_matrix./max(radar_matrix(:)))));
    figure
    imagesc(Range,time_axis,RTI_plot,[-50,0]); 
    xlim([0 100]);
    grid on            
    colorbar
    ylabel('Time (Sec)')
    xlabel('Range')      
    title('Monostatic RTI');

