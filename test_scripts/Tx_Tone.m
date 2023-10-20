clear all
close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
capture_duration = 6000;    % capture duration
Fs = 30e6;               % Sample Rate of SDR per I & Q (in reality Fs is double this)
pulse_duration = 2e-3;   % Desired Pulse Duration 
Bw = 30e6;               % LFM Bandwidth 

% Radar Parameters 
Fc = 2.44e9;   % Central RF 
Tx_gain = 66;  % [-23.75, 66] (S-Band = 23.5 dBm) (C-Band = 15.8 dBm)
Tx_SDR = 10;   % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel

% Parameters not configurable by user 
    C = physconst('LightSpeed');
    PRF = 1/pulse_duration;
    slope = Bw/pulse_duration;
    F_Max = Bw/2;
    R_Max = beat2range(F_Max,slope);
    sample_duration = 1/Fs;
    samples_per_pulse = pulse_duration/sample_duration;
    number_pulses = capture_duration / pulse_duration;
    number_cap_samps = capture_duration/sample_duration;
    Fc_M = Fc/1e6;   % RF in MHz 
    Bw_M = Bw/1e6;      % BW in MHz


%% Create Sawtooth Chirp for bladeRF
chirp = tone_generator(14.5e6,pulse_duration,Fs);
% figure; plot(real(chirp))

chirp = chirp;%.*window('blackman',size(chirp,1));
save_sc16q11('/tmp/chirp.sc16q11', chirp); %save chirp to binary file
% clear chirp
%     figure
%     spectrogram(chirp,128,100,128,Fs,'centered','yaxis') %plot spectrogram of chirp 
%     f = linspace(-0.5 * Fs, 0.5 * Fs, length(chirp));
%     figure
%     plot(f,20*log10(abs(fftshift(fft(chirp)))/samples_per_pulse));%plot FFT of chirp to show entire spectral content

%% Setup SDR
    % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''

    % Setup Tx SDR 
    [trig_flag_1,tx_command] = create_shell_command( ...
                                   false,...
                                   Experiment_ID,...
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
                                   3,...
                                   'tx');
    %tx_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(tx_command);

  