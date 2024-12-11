clear all
addpath('~/Documents/bladeRAD-Rikki/generic_scripts/matlab',...
        '~/Documents/bladeRAD-Rikki/generic_scripts',...
        '~/Documents/bladeRAD-Rikki/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 0005;       % Expeiment Name
capture_duration = 2;    % capture duration
Bw = 25e6;               % Sample Rate of SDR per I & Q (in reality Fs is double this)30 no overun 40 has overun
save_directory = "~/Documents/bladerad_data/Captures/"; % each experiment will save as a new folder in this directory
passive_max_range = 1000; %max range to cross-correlate tco

% Radar Parameters 
Fc = 2225e6;   % Central RF    
Ref_gain = 47;
Sur_gain = 60;
Pass_SDR = 1;   % SDR to use for Passive Radar - labelled on RFIC Cover and bladeRAD Facia Panel

% Parameters not configurable by user 
    C = physconst('LightSpeed');
    Fs = Bw;
    sample_duration = 1/Fs;
    number_cap_samps = 2*(capture_duration/sample_duration);
    RF_freq = Fc/1e6;   % RF in MHz 
    Bw_M = Bw/1e6;      % BW in MHz
    file_size_MBytes = (number_cap_samps * 16)*2/(8*1e6) ;

    
%% Setup Radar
    % 1 'set clock_sel external'; 2 'set clock_out enable; 3 'set clock_ref enable'

    % Setup Passive SDR 
   [trig_flag, passive_command] = create_shell_command(true,...
                                   Experiment_ID,...
                                   number_cap_samps,... 
                                   0,...
                                   0,...
                                   0,...
                                   Ref_gain,...
                                   Sur_gain,...
                                   RF_freq,...
                                   Bw_M,...
                                   Pass_SDR,...
                                   'master',...
                                   3,...
                                   'rx');
                               
    %passive_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(passive_command);


    %% Save Raw Data and create header file for directory 
    exp_dir = save_directory + Experiment_ID + '/';
    make_dir = 'mkdir ' + exp_dir;
    system(make_dir); % Blocking system command execution
    move_file = 'mv /tmp/active_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Rx Data Copyied to Save directory"
    else 
        "Rx Copy Failed"
        return
    end
    save(exp_dir + 'Passive Experimental Configuration') 
