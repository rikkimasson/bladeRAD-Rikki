clear all
% close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...close
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 1088;    % Expeiment Name
capture_duration = 1;        % capture duration
save_directory = "~/Documents/bladeRAD_Captures/lab/"; % each experiment will save as a new folder in this directory
exp_dir = save_directory + Experiment_ID + '/';


% FMCW Parameters 
active.Fs = 30e6;          % Sample Rate of SDR per I & Q (in reality Fs is double this)
active.pulse_duration = 0.2e-3;   % Desired Pulse Duration 
active.Bw = 30e6;          % LFM Bandwidth 
active.Fc = 690e9;   % Central RF 
active.Tx_gain = 50;       % [-23.75, 66] (S-Band = 23.5 dBm) (C-Band = 15.8 dBm)
active.Rx1_gain = 60;      % [-16, 60]
active.Rx2_gain = 60;       % [-16, 60]
Rx_1_lna = true;
Rx_2_lna = true;
active.Tx_SDR = 1;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
active.Rx_SDR = 2;   % SDR to use for RX
% Procesing Parameters
active.max_range = 1.4e5; %max range to LPF filter data to

% Passive Radar Parameters
passive.Fc = 690e6;   % Central RF    sample_duration
passive.Ref_gain = -15; % 26dB seems good for C-Band Patch Antennas [-16, 60]
passive.Sur_gain = -15; % 40dB seems good for C-Band Patch Antennas [-16, 60]
passive.SDR = 3;   % SDR to use for Passive Radar - labelled on RFIC Cover and bladeRAD Facia Panel
passive.Bw = 10e6;
passive.Fs = passive.Bw;
passive.max_range = 1000; %max range to cross-correlate to

process_active = true;
process_passive = false;
 
% Parameters not configurable by user 
    C = physconst('LightSpeed');
    % FMCW parameters
    active.PRF = 1/active.pulse_duration;
    active.slope = active.Bw/active.pulse_duration;
    active.F_Max = active.Bw/2;
    R_Max = beat2range(active.F_Max,active.slope);
    active.sample_duration = 1/active.Fs;
    active.active.samples_per_pulse = active.pulse_duration/active.sample_duration;
    active.number_pulses = capture_duration / active.pulse_duration;
    active.number_cap_samps = capture_duration/active.sample_duration;
    active.RF_freq = active.Fc/1e6;   % RF in MHz 
    active.Bw_M = active.Bw/1e6;      % BW in MHz
    
    % passive parameters
    passive_Fs = passive.Bw;
    passive_sample_duration = 1/passive_Fs;
    passive_number_cap_samps = (2*capture_duration)/passive_sample_duration;
    passive_RF_freq = passive.Fc/1e6;   % RF in MHz 
    passive.Bw_M = passive.Bw/1e6;      % BW in MHz



% %% Setup Passive Radar
    % by default clock out is always enabled 
    % 1 'set clock_sel external'; 
    % 2 'set clock_ref enable';

    % Setup Passive SDR 
   [trig_flag_3, passive_command] = create_shell_command(false,...
                                   Experiment_ID,...
                                   passive_number_cap_samps,... 
                                   0,...
                                   0,...
                                   0,...
                                   passive.Ref_gain,...
                                   passive.Sur_gain,...
                                   passive_RF_freq,...
                                   passive.Bw_M,...
                                   passive.SDR,...
                                   'master',...
                                   3,...
                                   'pass');
                               

    %passive_command = passive_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(passive_command);





    
    
    
    
 %% Save Raw Data and create header to directory 
    % make save directory
        make_dir = 'mkdir ' + exp_dir;
        system(make_dir);
    % move passive file to save directory
    move_file = 'mv /tmp/passive_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Passive Data Copied to Save directory"
    else 
        "Passive Copy Failed"
        return
    end
    save(exp_dir + 'Experimental Configuration') 


 %% Passive Processing
                            % survallance and reference. Will affect SNR dramatically.
    % load signal and split ref and sur
        file_location = exp_dir + 'passive_' + Experiment_ID;
        [ref_channel, sur_channel]  = load_passive_data(file_location);
    % Plot time domain signals
%          figure
%          fig = subplot(2,1,1);
%             plot(real(ref_channel(1:4000000)));
%             title("Ref channel time series");
%          hold on
%          subplot(2,1,2)
%              plot(real(sur_channel(1:4000000)));
%              title("Sur channel time series");    
%              fig_name = exp_dir + "Time Domain Signals_" + Experiment_ID + ".jpg";
%              saveas(fig,fig_name,'jpeg');


             nfft= 1000;
% figure
% f = linspace(-0.5 * passive.Fs, 0.5 * passive.Fs, nfft);
% freq_dom = fftshift(abs(fft(ref_channel(1:nfft),nfft)));
% % plot(f,20*log10(freq_dom/sqrt(50))-pow2db(nfft))
% plot(f,20*log10(freq_dom)- 20*log10(nfft-2^12/2)+10)
% grid on; grid minor;
% xlabel('Frequency');
% ylabel('Relative Power')

adc_max = max(real(sur_channel))
adc_rms = rms(real(sur_channel))
dbfs = 10*log10(adc_max/1)
dbm = 30+ 20*log10(adc_rms/(sqrt(50)))

raw_data = (ref_channel);

nfft= 4096;
f = linspace(-0.5 * passive.Fs, 0.5 * passive.Fs, nfft);
freq_dom = abs(pwelch(raw_data(1:nfft),nfft,0.8,'centered','power'));
pwr_lna_RX_2 = 30+10*log10(freq_dom)- 10*log10(nfft+ 2^12/2)+19;
plot(f,pwr_lna_RX_2); hold on;
plot(f,pwr_lna_RX_1);

grid on; grid minor;
xlabel('Frequency');
ylabel('Relative Power')

dbfs = 10*log10(adc_rms*sqrt(2)/1)
dbm = 30+ 20*log10(adc_rms/(sqrt(50)))

peak_freq=max(30+10*log10(freq_dom)- 10*log10(nfft+ 2^12/2)+19)