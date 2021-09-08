addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 666;    % Expeiment Name
Cap_dur = 4;        % capture duration
Fs = 40e6;          % Sample Rate of SDR per I & Q (in reality Fs is double this)
Pulse_dur = 1e-3;   % Desired Pulse Duration 
Bw = 40e6;          % LFM Bandwidth 

% Radar Parameters 
Fc = 500e6;         % Central RF 
Tx_gain = 10;       
Rx_gain = 36;
TX_SDR = 1;         % SDR to use for TX
RX_SDR = 2;         % SDR to use for RX

% Parameters not configurable by user 
    C = physconst('LightSpeed');
    PRF = 1/Pulse_dur
    Slope = Bw/Pulse_dur;
    F_Max = Bw/2;
    R_Max = beat2range(F_Max,Slope);
    Samp_dur = 1/Fs;
    Samps_per_pulse = Pulse_dur/Samp_dur;
    Num_pulses = Cap_dur / Pulse_dur;
    Num_cap_samps = Cap_dur/Samp_dur;
    RF_freq = Fc/1e6;   % RF in MHz 
    Bw_M = Bw/1e6;      % BW in MHz


%% Create Sawtooth Chirp for bladeRF
chirp = saw_LFM_chirp(Bw,Pulse_dur,Fs);
save_sc16q11('/tmp/chirp.sc16q11', chirp); %save chirp to binary file
    % spectrogram(chirp,128,100,128,Fs,'centered','yaxis') %plot spectrogram of chirp
    % f = linspace(-0.5 * Fs, 0.5 * Fs, length(chirp));
    % figure
    % plot(f,20*log10(abs(fftshift(fft(chirp)))/Samps_per_pulse));%plot FFT of chirp to show entire spectral content

%% Setup Radar
% 1 'set clock_sel external' 
% 2 'set clock_out enable'
% 3 'set clock_ref enable'

% Setup Tx SDR - waiting for external trigger from RX SDR / output clk
tx_command = create_shell_command(Experiment_ID,...
                               Num_cap_samps,... 
                               Num_pulses,...
                               Tx_gain,...
                               Rx_gain,...
                               RF_freq,...
                               Bw_M,...
                               TX_SDR,...
                               'slave',...
                               2,...
                               'tx');

%tx_command = '/home/piers/repos/bladeRAD/generic_scripts/' + tx_command + "&"; % Non-Blocking system command execution                    
status = system(tx_command);

% % Setup Rx SDR - trigger Tx / disapline from Tx SDR Clock 
% rx_command = create_shell_command(Experiment_ID,...
%                                Num_cap_samps,... 
%                                Num_pulses,...
%                                Tx_gain,...
%                                Rx_gain,...
%                                RF_freq,...
%                                Bw_M,...
%                                RX_SDR,...
%                                'master',...
%                                'set clock_sel external',...
%                                'rx');
%         
% system(rx_command); % Blocking system command execution 







% 
% % Load Reference Deramp Signal
% refsig = load_refsig(Bw,Fc,Pulse_dur);
% % figure 
% % spectrogram(refsig,128,100,100,Fs,'centered','yaxis')
% 
% 
% %% Load Signal    
% %%load file in 30s chunks
% window_size = 30; %window size in s
% windows = Cap_dur / window_size;
% 
% for window=0:windows-1
%     
%     %load rx file
%     tic
% %     rx_file= "/home/piers/Documents/bladeRF_Code/FMCW/Experiments/" + Test_id + ".sc16q11";
%     rx_file= "/media/piers/data_drive/BladeRF_Experiments/Captures/" + Test_id + "/N1_" + Test_id + ".sc16q11"
% %     raw_data= load_sc16q11(rx_file); 
%     raw_data = load_sc16q11_window(rx_file,window,Fs,window_size);
%     loading_rx_file = toc 
% % plot(real(raw_data)) check the bits toggled at ADC
% 
% %% Reshape array into matrix of pulses
%     tic
%     pulse_matrix = reshape(raw_data,[length(raw_data)/(Num_pulses/windows),Num_pulses/windows]); %reshape array to individual pulses
%     clear raw_data
% %     figure 
% %     spectrogram(pulse_matrix(:,1000),128,100,100,Fs,'centered','yaxis')
%     reshape_radar_matrix = toc 
% 
% 
% %% Range Limmited Signal Processing 
% 
% 
% %% Deramp Signal
%     Data_Deramped = zeros(size(pulse_matrix,1),size(pulse_matrix,2));
%     tic
%     Data_Deramped = dechirp(pulse_matrix,refsig);
%     clear pulse_matrix
%     deramp_data = toc
% 
% 
% 
% %% Decimate Signal 
% %%determine decimation factor 
% %decimation factor must be by an integer thus the maximum range will likley
% %be non integer due to rounding
%     d_max = 1000; % set an approximate maximum range in metres
%     if_max = (d_max * 2 * Slope)/C; % if_max is the maximum frequency were intrested in
%     dec_fac = Fs/if_max; % dec_fac is the non integer ideal decimation factor
%     dec_fac_rnd = round(size(Data_Deramped,1)/dec_fac); %dec_fac_rnd is the closest integer 
%     dec_fac_actual = round(size(Data_Deramped,1)/dec_fac_rnd); %determines actual optmal decimation factor
%     if_freq_actual = Fs/dec_fac_actual; %determines subsequent maximum if frequency 
%     d_max_actual = beat2range(if_freq_actual,Slope); %translates if frequency to range
% 
% 
% %decimate signal
%     Dec_Deramped = zeros(dec_fac_rnd + 1,size(Data_Deramped,2)); %initiate array
%     tic
%     for i=1:Num_pulses/windows
%     Dec_Deramped(:,i) = decimate(Data_Deramped(:,i),dec_fac_actual);
%     end
%     clear Data_Deramped
%     Decimation_time = toc
% 
% %% FFT Signal
%     tic
%     %         Dec_Deramped = fft(Dec_Deramped,blackman(size(Dec_Deramped,1)),size(Dec_Deramped,1)); %Convert to Power Spectral Density
% %     Final_Data = fft(Dec_Deramped,size(Dec_Deramped,1));
%     windowed_data{window+1} = fft(Dec_Deramped,size(Dec_Deramped,1));
%     clear Dec_Deramped 
% end
% 
% %% Save data to file
% save_directory = "/media/piers/data_drive/BladeRF_Experiments/Captures/" + Test_id;
% sys_command = "mkdir " + save_directory;
% % system(sys_command)
% mat_name = save_directory + "/N1_" + Test_id;
% save(mat_name,'windowed_data','Num_pulses','Fs','Pulse_dur')
% 
% Final_Data  = [windowed_data{1}, windowed_data{2}, windowed_data{3}, windowed_data{4}];
% 
% save_directory = "/media/piers/data_drive/BladeRF_Experiments/Captures/" + Test_id;
% sys_command = "mkdir " + save_directory;
% % system(sys_command)
% mat_name = save_directory + "/Final_Data_Backup" + Test_id;
% save(mat_name,'Final_Data','Num_pulses','Fs','Pulse_dur')
% 
% % 
% % Bi_Matched_Data = Final_Data;
% 
% % Mono_Final_Data = Final_Data;
% 
% % plot(real(raw_data));
% 
% %% Plot RTI
% 
% tic 
% % Range_axis = linspace(0,d_max_actual,size(Final_Data,1));
% Range_bin = 1:size(Final_Data,1);
% time_axis = linspace(0,size(Final_Data,2)*Pulse_dur,size(Final_Data,2));
% RTI_plot= transpose(10*log10(abs(Final_Data./max(Final_Data(:)))));
% 
% % 
% % %Plot RTI over Range in meters 
% % figure
% % imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
% % % xlim([0 25])
% % %ylim([0 0.0005])
% % grid on            
% % colorbar
% % ylabel('Time (Sec)')
% % xlabel('Range (m)')   
% % fig_title = "Monostatic RTI - Test " + Test_id;
% % title(fig_title);fig_name = save_directory + "/RTI_" + Test_id + ".jpg";
% % saveas(fig,fig_name,'jpeg')
% 
% % % fig_name = "/home/piers/Desktop/FMCW/Experiments/Test_"+ Test_id + ".fig";
% % % savefig(fig_name)
% % plot_signal = toc
% 
% %Plot RTI over Range in Range Bins
% figure
% fig = imagesc(Range_bin,time_axis,RTI_plot,[-50,0]);
% xlim([1 20])
% %ylim([0 0.0005])
% grid on            
% colorbar
% ylabel('Time (Sec)')
% xlabel('Range Bin')   
% fig_title = "Monostatic RTI - Test " + Test_id;
% title(fig_title);
% fig_name = save_directory + "/RTI_" + Test_id + ".jpg";
% saveas(fig,fig_name,'jpeg')
% plot_signal = toc
% % 
% % %% Coherent integration 
% % compressed_data = sum(Dec_Deramped,2);
% % figure
% % plot(Range_axis,abs(compressed_data))
% % xlim([0 1000])
% 
% % %% Spectrogram 
% % spec_bin = 3;
% % l_fft = 512;
% % pad_factor = 4;
% % overlap_factor = 0.99;
% % [spect,f] = spectrogram(Final_Data(spec_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
% % % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
% % v=dop2speed(f,C/Fc)*2.237;
% % spect= 10*log10(abs(spect./max(spect(:))));
% % figure
% % fig = imagesc(time_axis,f,spect,[-50 0]);   
% % ylim([-100 100])
% % colorbar
% % xlabel('Time (Sec)')
% % % ylabel('Radial Velocity (mph)')   
% % ylabel('Doppler Frequency (Hz)')  
% % fig_title = "Monostatic Spectrogram - Test " + Test_id;
% % title(fig_title);
% % fig_name = save_directory + "/Spectrogram_" + Test_id + ".jpg";
% % saveas(fig,fig_name,'jpeg')
% % 
% % 
% 
% 
% 
% phase = unwrap(angle(Final_Data(6,:)));
% [psd_Mono,f] = pwelch(phase,Num_pulses,0,Num_pulses,PRF,'centered','power');
% % [psd_Mono,f] = pwelch(phase,Num_pulses,512,2048,PRF,'centered','power');
% phase = phase - max(phase(100));
% phase = rad2deg(phase);
% phase = (phase/360) * (1/2.4e9);
% 
% Bi_phase = unwrap(angle(Bi_Matched_Data(6,:)));
% [psd_Bi,f] = pwelch(Bi_phase,Num_pulses,0,Num_pulses,PRF,'centered','power');
% % [psd_Bi,f] = pwelch(Bi_phase,Num_pulses,512,2048,PRF,'centered','power');
% Bi_phase = Bi_phase - max(Bi_phase(100));
% Bi_phase = rad2deg(Bi_phase);
% Bi_phase = (Bi_phase/360) * (1/2.4e9) +1.5e-11;
% 
% 
% figure
% plot(time_axis,phase)
% hold on
% plot(time_axis,Bi_phase)
% grid on            
% ylabel('Phase (s)')
% xlabel('Time (s)')   
% fig_title = "Monostatic Phase of 2.4GHz Pulse";
% legend('Common VCTCXO','Common 10MHz')
% % title(fig_title);
% 
% figure
% p1 = semilogx(f,10*log10(psd_Mono));
% hold on
% p2 = semilogx(f,10*log10(psd_Bi));
% p1.Color(4) =0.5;
% p2.Color(4) =0.2;
% xlim([0 2500]);
% ylabel('dBrads^{2}/Hz')
% xlabel('Frequency (Hz)')   
% legend('Common VCTCXO','Common 10MHz')
% grid on     
% %     spectrogram(refsig,128,100,100,Fs,'centered','yaxis')
% 
% doppler_info = fftshift(fft(Final_Data(7,1:end-1),Num_pulses-1));
% f = linspace(-PRF/2,PRF/2,Num_pulses-1);
% figure
% plot(f,10*log10(abs(doppler_info)));
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
% 
% 
% 
% 
% 
% 
% % %% MTI Filtering 
% % % Single Delay Line Filter 
% % MTI_Data = zeros(size(Final_Data));
% %       for i=2:Num_pulses
% %             MTI_Data(:,i) = Final_Data(:,i)-Final_Data(:,i-1);
% %       end
% %       
% % %Plot MTI RTI      
% % MTI_RTI_plot= transpose(10*log10(abs(MTI_Data./max(MTI_Data(:)))));
% % figure
% % fig = imagesc(Range_bin,time_axis,MTI_RTI_plot,[-50,0]);
% % xlim([1 20])
% % %ylim([0 0.0005])
% % grid on            
% % colorbar
% % ylabel('Time (Sec)')
% % xlabel('Range Bin')   
% % fig_title = "Monostatic Single Delay Line MTI  RTI - Test " + Test_id;
% % title(fig_title);
% % fig_name = save_directory + "/MTI_RTI_" + Test_id + ".jpg";
% % saveas(fig,fig_name,'jpeg')
% % plot_signal = toc     
% % 
% % 
% % %Plot MTI Spectrogram  
% % [spect,f] = spectrogram(MTI_Data(spec_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
% % % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
% % v=dop2speed(f,C/Fc)*2.237;
% % spect= 10*log10(abs(spect./max(spect(:))));
% % figure
% % fig = imagesc(time_axis,f,spect,[-50 0]);
% % ylim([-100 100])
% % colorbar
% % xlabel('Time (Sec)')
% % % ylabel('Radial Velocity (mph)')   
% % ylabel('Doppler Frequency (Hz)')  
% % fig_title = "Monostatic Single Delay Line MTI Spectrogram - Test " + Test_id;
% % title(fig_title);
% % fig_name = save_directory + "/MTI_Spectrogram_" + Test_id + ".jpg";
% % saveas(fig,fig_name,'jpeg')
% 
