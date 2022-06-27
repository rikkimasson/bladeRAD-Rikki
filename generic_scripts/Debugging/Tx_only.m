clear all
close all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 8;       % Expeiment Name
capture_duration = 30  ;    % capture duration
Fs = 20e6;               % Sample Rate of SDR per I & Q (in reality Fs is double this)
pulse_duration = 1e-3;   % Desired Pulse Duration 
Bw = 20e6;               % LFM Bandwidth 
save_directory = "/home/piers/Documents/Captures/7_Oct/FMCW/"; % each experiment will save as a new folder in this directory
    exp_dir = save_directory + Experiment_ID + '/';

% Radar Parameters 
Fc = 2.4e9;   % Central RF 
Tx_gain = 0;  % [-23.75, 66] (S-Band = 23.5 dBm) (C-Band = 15.8 dBm)
Rx1_gain = 16;  % [-16, 60]
Rx2_gain = 0;  % [-16, 60]
Tx_SDR = 3;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
Rx_SDR = 2;   % SDR to use for RX

% Procesing Parameters
 max_range = 1000; %max range to LPF filter data to

% Parameters not configurable by user 
    C = physconst('LightSpeed');
    PRF = 1/pulse_duration;
    slope = Bw/pulse_duration;
    F_Max = Bw/2;
    R_Max = beat2range(F_Max,slope);
    sample_duration = 1/Fs;
    samples_per_pulse = pulse_duration/sample_duration
    number_pulses = capture_duration / pulse_duration
    number_cap_samps = capture_duration/sample_duration
    Fc_M = Fc/1e6;   % RF in MHz 
    Bw_M = Bw/1e6;      % BW in MHz


%% Create Sawtooth Chirp for bladeRF
chirp1 = saw_LFM_chirp(Bw,pulse_duration,Fs);
% figure
save_sc16q11('/tmp/chirp.sc16q11', chirp1); %save chirp to binary file
clear chirp
% spectrogram(chirp,128,100,128,Fs,'centered','yaxis') %plot spectrogram of chirp
% f = linspace(-0.5 * Fs, 0.5 * Fs, length(chirp));

%% Setup Radar
    % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''

    % Setup Tx SDR 
    [trig_flag_1,tx_command] = create_shell_command(Experiment_ID,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   0,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   Fc_M,...
                                   Bw_M,...
                                   Tx_SDR,...
                                   'master',...
                                   1,...
                                   'tx');
%     tx_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(tx_command);
   
    
    pause(5);
% 
%     % Setup Rx SDR 
%     [trig_flag_2,rx_command] = create_shell_command(Experiment_ID,...
%                                    number_cap_samps,... 
%                                    number_pulses,...
%                                    0,...                                 
%                                    Tx_gain,...
%                                    Rx1_gain,...
%                                    Rx2_gain,...
%                                    Fc_M,...
%                                    Bw_M,...
%                                    Rx_SDR,...
%                                    'master',...
%                                    1,...
%                                    'rx'); 
%     if trig_flag_1 && trig_flag_2
%         "Trigger Conflict - FMCW Radar"
%         return
%     end                                                          
%     system(rx_command) % Blocking system command execution 


