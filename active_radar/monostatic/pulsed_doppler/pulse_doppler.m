%% Capture Parameters
C = physconst('LightSpeed');
Fs = 40e6;          %Sample Rate of SDR per I & Q (in reality Fs is double this)
Pulse_dur = 100e-6;   %Desired Pulse Duration 
Bw = 40e6;          %LFM Bandwidth 
Samp_dur = 1/Fs;
Samps_per_pulse = Pulse_dur/Samp_dur;

Num_pulses = 100
PRF = 1000

PRI = 1/PRF;
Samps_per_PRI = PRI/Samp_dur;
Cap_dur = Num_pulses / PRF
Num_cap_samps = Cap_dur/Samp_dur;
Tx_Delay = PRI - Pulse_dur % in seconds
Tx_Delay_us = Tx_Delay * 10e5 % in micro seconds
R_Max = PRI*C/2;

%% Create Sawtooth Chirp for bladeRF
  chirp = SAW_LFM_CHIRP(Bw,Pulse_dur,Fs);
% spectrogram(chirp,128,100,128,Fs,'centered','yaxis') %plot spectrogram of chirp
% f = linspace(-0.5 * Fs, 0.5 * Fs, length(chirp));
% figure
% plot(f,20*log10(abs(fftshift(fft(chirp)))/Samps_per_pulse));%plot FFT of chirp to show entire spectral content
save_sc16q11('/tmp/40MHz_chirp.sc16q11', chirp); %save chirp to binary file

%% Radar Parameters 

Tx_gain = 35;
Rx_gain = 0;
Test_id = 6;
RF_freq = 500; %in MHz 
Bw_M = Bw/1e6;

%Transeiver Command Generation
N0_command = "/home/piers/Documents/bladeRF_Code/PulseDoppler/Shell_Scripts/bi_N0_PD.sh " + Num_cap_samps + " " + Tx_Delay_us + " " + Num_pulses + " " + Tx_gain + " " + Rx_gain + " " + RF_freq + " " + Bw_M + " " + Test_id; 
%Passive Receiver Command Generation 
N1_command = "/home/piers/Documents/bladeRF_Code/PulseDoppler/Shell_Scripts/N1_PD.sh " + Num_cap_samps + " " + Rx_gain + " " + RF_freq + " " + Bw_M + " " + Test_id + " &"; 

status = system(N1_command);
pause(5);
status = system(N0_command);


%% Load Signals
N0_rx_file = "/home/piers/Documents/bladeRF_Code/PulseDoppler/Experiments/N0_" + Test_id + ".sc16q11";
N1_rx_file = "/home/piers/Documents/bladeRF_Code/PulseDoppler/Experiments/N1_" + Test_id + ".sc16q11";
%rx_file= "/Users/piersbeasley/OneDrive - University College London/BladeRF System/FMCW/Experiments/" + Test_id + ".sc16q11";
N0_raw_data = load_sc16q11(N0_rx_file); 
N1_raw_data = load_sc16q11(N1_rx_file); 

%% Reshape array into matrix of pulses
N0_pulse_matrix = reshape(N0_raw_data,[length(N0_raw_data)/Num_pulses,Num_pulses]); %reshape array to individual pulses
N1_pulse_matrix = reshape(N1_raw_data,[length(N1_raw_data)/Num_pulses,Num_pulses]); %reshape array to individual pulses

spectrogram(N0_raw_data(1:4000),128,100,100,Fs,'centered','yaxis')


% Load Reference Deramp Signal
if RF_freq == 500
   refsig = load_sc16q11("/home/piers/Documents/bladeRF_Code/PulseDoppler/MATLAB/ref_waveforms/100us_500MHz.sc16q11");
else 
    refsig = load_sc16q11("/home/piers/Documents/bladeRF_Code/PulseDoppler/MATLAB/ref_waveforms/1ms_2.4g.sc16q11");
end


%% Range Limmited Signal Processing 

% Match Filter 
% refsig = transpose(refsig);
refsig = refsig(1:Samps_per_pulse);


refsig_fft = conj(fft(refsig,Samps_per_PRI)); 
N0_Data_matched = zeros(size(N0_pulse_matrix));

for i = 1:Num_pulses
appo = fft(N0_pulse_matrix(1:end,i));
N0_Data_matched(:,i) = ifft(appo.*refsig_fft); 
end

%Range_bin = 1:size(Data_matched,1) * samp_dur * C;
Range_bin = 1:size(N0_Data_matched,1);
Range = 1:R_Max;
time_axis = 1:Cap_dur;%linspace(0,PRI,size(N0_Data_matched,2));
RTI_plot= transpose(10*log10(abs(N0_Data_matched./max(N0_Data_matched(:)))));
figure
imagesc(Range,time_axis,RTI_plot,[-50,0]); 
 xlim([0 300]);
grid on            
colorbar
ylabel('Time (Sec)')
xlabel('Range')      
title('Monostatic RTI');



%% bistatic capture


N1_Data_matched = zeros(size(N1_pulse_matrix));

for i = 1:Num_pulses
appo = fft(N1_pulse_matrix(1:end,i));
N1_Data_matched(:,i) = ifft(appo.*refsig_fft); 
end

Range_bin = 1:size(N1_Data_matched,1);
Range = 1:R_Max;
time_axis = 1:Cap_dur;%linspace(0,PRI,size(N1_Data_matched,2));
RTI_plot= transpose(10*log10(abs(N1_Data_matched./max(N1_Data_matched(:)))));
figure
imagesc(Range_bin,time_axis,RTI_plot,[-50,0]); 
xlim([0 300]);
grid on            
colorbar
ylabel('Time (Sec)')
xlabel('Range')      
title('Bistatic RTI');

f = linspace(-0.5 * PRF, 0.5 * PRF ,size(N1_Data_matched,2)-1);
figure

sum_rbins = sum(N1_Data_matched(40:55,2:end));

plot(f,20*log10(abs(fftshift(fft(sum_rbins)))));
ylabel('Relative Power (dB)')
xlabel('Frequency Offset (Hz)')  
title('RF Frequency Offset Between Monostatic and Bistatic Node');

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % plot(real_rx)
% 
% % f = linspace(-0.5 * Fs, 0.5 * Fs, length(rx_signal));
% % figure
% % plot(f,20*log10(abs(fftshift(fft(rx_signal)))/NUM_CAP_SAMPLES));
% % xlabel('Frequency (Hz)');
% % ylabel('Power (dB)');
% % title('Bi_20 MHz LFM; 20MHz Fs');
% % figure
% % spectrogram(rx_signal,128,100,100,Fs,'c-entered','yaxis')
% % %xlim([50 60]);
% % title('Bi_15 MHz LFM');
% 
% % Match filter results
% % reshape array into matrix of pulses
% raw_data = bi_rx_signal;
% pulse_matrix = reshape(raw_data,[SAMPS_PER_PRI,NUM_PULSES]); %reshape array to individual pulses
% figure
% spectrogram(raw_data(1: 2000000),128,100,100,Fs,'centered','yaxis')
% %spectrogram(pulse_matrix(:,15),128,100,100,Fs,'centered','yaxis')
% 
% % Match Filter the  
% ref_sig = transpose(refsig);% load_sc16q11('mon_Rx_20MHz_chirp_1e-4.sc16q11');
% ref_sig = Ref_sig(1:Samps_per_pulse);
% % figure
% %spectrogram(Ref_sig,128,100,100,Fs,'centered','yaxis')
% refsig_fft = conj(fft(Ref_sig,SAMPS_PER_PRI)); 
% Data_matched = zeros(size(pulse_matrix));
% 
% for i = 1:NUM_PULSES 
% appo = fft(pulse_matrix(1:end,i));
% Data_matched(:,i) = ifft(appo.*refsig_fft); 
% end
% %Range_bin = 1:size(Data_matched,1) * samp_dur * C;
% Range_bin = 1:size(Data_matched,1);
% time_axis = linspace(0,PRI,size(Data_matched,2));
% RTI_plot= transpose(10*log10(abs(Data_matched./max(Data_matched(:)))));
% figure
% imagesc(Range_bin,time_axis,RTI_plot,[-50,0]);
% %xlim([8070 8100]);
% grid on            
% colorbar
% ylabel('Time (Sec)')
% xlabel('Range Bin')      
% title('Bistatic RTI');
% f = linspace(-0.5 * PRF, 0.5 * PRF, size(Data_matched,2)-1);
% figure
% plot(f,20*log10(abs(fftshift(fft(Data_matched(8084,2:end))))/size(Data_matched,2)-1));
% ylabel('Relative Power (dB)')
% xlabel('Frequency Offset (Hz)')  
% title('RF Frequency Offset Between Monostatic and Bistatic Node');