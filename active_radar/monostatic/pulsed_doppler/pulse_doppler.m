clear all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 66;    % Expeiment Name
capture_duration = 10;        % capture duration
save_directory = "~/Documents/bladeRAD_Captures/15_09_2022_farm/pulse_doppler/"; % rach experiment will save as a new folder in this directory
exp_dir = save_directory + Experiment_ID + '/';

% Radar Parameters
Fc = 2.44e9;
Fs = 30e6;               %Sample Rate of SDR per I & Q (in reality Fs is double thiclos)
pulse_duration =  0.2e-3;  %Desired Pulse Duration 
Bw = 30e6;               %LFM Bandwidth 
PRF = 500;
Tx_gain = 66;   % [-23.75, 66] (S-Band = 23.5 dBm) (C-Band = 15.8 dBm)  
Rx1_gain = 18;   % [-16, 60]
Rx2_gain = 0;   % [-16, 60] 
Tx_SDR = 1;     % SDR to use for TX - labelled on RFIC Cover and bladeRAD Facia Panel
Rx_SDR = 2;     % SDR to use for RX

% Parameters not configurable by user 
    C = physconst('LightSpeed');
    PRI = 1/PRF;
    if PRI < pulse_duration
        "pulse duration longer than PRI"
        return
    end
    number_pulses = capture_duration/PRI;
    sample_duration = 1/Fs;
    samples_per_pulse = pulse_duration/sample_duration;
    samples_per_PRI = PRI/sample_duration;
    number_cap_samps = capture_duration/sample_duration
    Tx_Delay = PRI - pulse_duration; % in seconds
    Tx_Delay_us = Tx_Delay * 1e6; % in micro seconds
    R_Max = PRI*C/2;
    Fc_M = Fc/1e6;      % RF in MHz 
    Bw_M = Bw/1e6;      % BW in MHz

%% Create Sawtooth Chirp for bladeRF
 tx_chirp = saw_LFM_chirp(Bw,pulse_duration,Fs);
 save_sc16q11('/tmp/chirp.sc16q11', tx_chirp); %save chirp to binary file
 clear chirp


%% Setup Radar
    % 1 'set clock_sel external'; 2 'set clock_ref enable; 3 ''
 
    % Setup Tx SDR 
    [trig_flag_1,tx_command] = create_shell_command(true,...
                                   Experiment_ID,...
                                   number_cap_samps,... 
                                   number_pulses,...
                                   Tx_Delay_us,...
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
    status = system(tx_command);
    pause(5);   

    % Setup Rx SDR 
    [trig_flag_2,rx_command] = create_shell_command(false,...
                                   Experiment_ID,...
                                   number_cap_samps,... 
                                   0,...
                                   0,...
                                   Tx_gain,...
                                   Rx1_gain,...
                                   Rx2_gain,...
                                   Fc_M,...
                                   Bw_M,...
                                   Rx_SDR,...
                                   'master',...
                                   1,...
                                   'rx'); 
    if trig_flag_1 && trig_flag_2
        "Trigger Conflict"
        return
    end
%     rx_command = rx_command + "&"; % uncomment for non-blocking system command execution                                                              
    system(rx_command); % Blocking system command execution 

 
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

  %% Load Reference Deramp Signal
    refsig = load_refsig(Bw_M,Fs,Fc,pulse_duration);
%     figure
%     plot(abs(real(refsig)))
%     title("Reference Signal")

%% Load Signals
    file_location = exp_dir + 'active_' + Experiment_ID + '.sc16q11';
    raw_data = load_sc16q11(file_location);

    
%% Reshape array into matrix of pulses
pulse_matrix = reshape(raw_data,[length(raw_data)/number_pulses,number_pulses]); %reshape array to individual pulses
% figure
% spectrogram(pulse_matrix(:,1),128,100,100,Fs,'centered','yaxis')
figure
fig = plot(abs(real(pulse_matrix(:,1))));
    ylabel('ADC Value (0-1)')
    xlabel('Samples')    
    title("Pulse Doppler - Pulse Time Series - " + Experiment_ID);
    fig_name = exp_dir + "Pulse_Time_Series_" + Experiment_ID + ".jpg";
    saveas(fig,fig_name,'jpeg')



%% Range Limited Signal Processing 
zero_padding = 3;

% Match Filter 
matched_filter = refsig;
matched_filter_fft = conj(fft(matched_filter,(samples_per_PRI*zero_padding))); 
processed_signal = zeros(size(pulse_matrix,2),size(matched_filter_fft,1));
for i = 1:number_pulses
    pri = fft(pulse_matrix(:,i),size(matched_filter_fft,1));
    processed_signal(i,:) = ifft(pri.*matched_filter_fft); 
end

%% PLot RTI
    Range_bin = 1:size(processed_signal,2);
    Range = linspace(0,R_Max*2,size(processed_signal,2));
    time_axis = linspace(0,capture_duration,size(processed_signal,1));
    RTI_plot=10*log10(abs(processed_signal./max(max(processed_signal))));
    figure
    fig = imagesc(Range,time_axis,RTI_plot,[-50,0]); 
        xlim([0 1000]);
        grid on            
        colorbar
        ylabel('Time (Sec)')
        xlabel('Range')      
        title("Pulse Doppler - Monostatic RTI - " + Experiment_ID);
        fig_name = exp_dir + "Monostatic_RTI_" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')
    
%% Plot 5 single pulses      
    fig = figure
    ranges_2_plot = floor(linspace(1,number_pulses,5));
    for i = ranges_2_plot
        plot(Range,RTI_plot(i,:));
        hold on
    end
        title("Single Pulses - " + Experiment_ID);
        xlim([0 200])
        grid on
        ylabel('Relative Power (dB)')
        xlabel('Range (m)')  
        legend("pulse " + ranges_2_plot(1), "pulse " + ranges_2_plot(2), "pulse " + ranges_2_plot(3),"pulse " + ranges_2_plot(4),"pulse " + ranges_2_plot(5))
        fig_name = exp_dir + "Single_Pulse" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg') 
      
%% MTI Filtering 
        % Single Delay Line Filter 
        MTI_Data = zeros(size(processed_signal));
              for i=2:number_pulses
                    MTI_Data(i,:) = processed_signal(i,:)-processed_signal(i-1,:);
              end        
        


  %% PLot RTI
    Range_bin = 1:size(MTI_Data,2);
    Range = linspace(0,R_Max*2,size(MTI_Data,2));
    time_axis = linspace(0,capture_duration,size(MTI_Data,1));
    RTI_plot=10*log10(abs(MTI_Data./max(max(MTI_Data))));
    figure
    fig = imagesc(Range,time_axis,RTI_plot,[-50,0]); 
        xlim([0 1000]);
        grid on            
        colorbar
        ylabel('Time (Sec)')
        xlabel('Range')      
        title("Pulse Doppler - Monostatic RTI - " + Experiment_ID);
        fig_name = exp_dir + "Monostatic_RTI_" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')
    
%% Spectrogram rom
        r_start = 1;
        r_stop = 10;
        l_fft = 512;
        pad_factor = 1;
        overlap_factor = 0.99;
        integrated_data = sum(processed_signal(:,r_start:r_stop),2);
        [spect,f] = spectrogram(integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
        v=dop2speed(f,C/Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-inf 0]);   
%             ylim([-100 100])
            colorbar
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "FMCW Spectrogram -" + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID;
            saveas(fig,fig_name)
            
%% MTI Spectrogram
            MTI_integrated_data = sum(MTI_Data(:,r_start:r_stop),2);
            [spect,f] = spectrogram(MTI_integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
            v=dop2speed(f,C/Fc)*2.237;
            spect= 10*log10(abs(spect./max(spect(:))));
            figure
            fig = imagesc(time_axis,f,spect,[-50 0]);   
                ylim([-100 100])
                colorbar
                xlabel('Time (Sec)')
                % ylabel('Radial Velocity (mph)')   
                ylabel('Doppler Frequency (Hz)')  
                fig_title = "MTI FMCW Spectrogram - " + Experiment_ID;
                title(fig_title);
                fig_name = exp_dir + "MTI FMCW Spectrogram_" + Experiment_ID + ".jpg";
                saveas(fig,fig_name,'jpeg')
                fig_name = exp_dir + "MTI FMCW Spectrogram_" + Experiment_ID;
                saveas(fig,fig_name)

