% clear all
%close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 500;       % Expeiment Name
capture_duration = 0.1;    % capture duration
Fs = 30e6;               % Sample Rate of SDR per I & Q (in reality Fs is double this)
pulse_duration = 1e-3;   % Desired Pulse Duration 
Bw = 30e6;               % LFM Bandwidth 
save_directory = "~/Documents/bladeRAD_Captures/lab/"; % each experiment will save as a new folder in this directory
exp_dir = save_directory + Experiment_ID + '/';

% Radar Parameters 
Fc = 2.44e9;   % Central RF 
Rx1_gain = 60;  % [-15, 60]
Rx2_gain = 0;  % [-6, 60]
Rx_SDR = 10;   % SDR to use for RX

% Procesing Parameters
 max_range = 1000; %max range to LPF filter data to
 processing_flag = true;

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


%% Setup SDR
    % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''

    % Setup Rx SDR 
    [trig_flag_2,rx_command] = create_shell_command( ...
                                   false,...
                                   Experiment_ID,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   0,...                                 
                                   0,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   Fc_M,...
                                   Bw_M,...
                                   Rx_SDR,...
                                   'master',...
                                   3,...
                                   'rx'); 
                                                       
    system(rx_command) % Blocking system command execution 


%% Save Raw Data and create  header to directory 
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
    save(exp_dir + 'Experimental Configuration') 

%% load rx file
file_location = exp_dir + 'active_' + Experiment_ID;
rx_file= file_location+ ".sc16q11";
raw_data = load_sc16q11(rx_file);

plot(real(raw_data)) %check the bits toggled at ADC
adc_max = max(real(raw_data))
adc_rms = rms(real(raw_data))
dbfs = 10*log10(adc_rms*sqrt(2)/1)
dbm = 30+ 20*log10(adc_rms/(sqrt(50)))


nfft= 4096;
f = linspace(-0.5 * Fs, 0.5 * Fs, nfft);
freq_dom = abs(pwelch(raw_data(1:nfft),nfft,0.8,'centered','power'));
pwr_lna_RX_2 = 30+10*log10(freq_dom)- 10*log10(nfft+ 2^12/2)+19;
plot(f,pwr_lna_RX_2);
grid on; grid minor;
xlabel('Frequency');
ylabel('Relative Power')



        