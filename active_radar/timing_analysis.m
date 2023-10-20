close all;
clear all;
% paths to generic functions
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts/matlab/CFAR/',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/') 

% Local directory containing the radar data
    local_save_directory =  "/media/sdrlaptop1/T7/07_08_2023_farm_n0/";

% Processing Flags 
    process_active_a = true;
    process_passive_a = false;
     
% Select the experiment you wish to process
    experiment_number = 1;

for i=experiment_number
    
%% Load .mat file containing experiment setup parameters
        mat_file_name = local_save_directory + i + "/Experimental Configuration.mat";
        load(mat_file_name);
        exp_dir = local_save_directory + i + '/';


%% FMCW Processing and Print RTI
if process_active_a == true
       % load refsig for deramping
        refsig = load_refsig(active.Bw_M,active.Fs,active.Fc,active.pulse_duration);    
    % load Signal, Mix and Dermap Signal  
        file_location = exp_dir + 'active_' + Experiment_ID;
        lp_filter = true;
%         active.max_range = 1000
        [max_range_actual,deramped_signal,active.decimation_factor_actual] = deramp_and_decimate(file_location,active.max_range,refsig,capture_duration,active.number_pulses,active.Fs,active.slope,lp_filter);
    
    % Window and FFT Signal 
    % window signal
        beat_frequncies = deramped_signal;%.*window('blackman',size(deramped_signal,1));
        clear  deramped_signal refsig

    % fft signal
        zero_padding = 1; % 1 = none; 2 = 100%
        beat_frequncies = fft(beat_frequncies,size(beat_frequncies,1)*zero_padding);
        beat_frequncies = beat_frequncies(1:(size(beat_frequncies,1)/2),:); % keep +ve beat frequencies

% Derive range and time axis 
        active.n_range_bins = size(beat_frequncies,1);
        active.range_bins = 1:active.n_range_bins;
        active.fftfrequncies =fftfreq(size(beat_frequncies,1),1/(active.Fs/active.decimation_factor_actual)); % possible beat frequencies
        active.slope = active.Bw/active.pulse_duration;
        ranges = (active.fftfrequncies*C)/(2*active.slope); % calculate true range bin size    
        active.range_axis = ranges(1:(size(ranges,2)/2)); % truncate to only +ve beat frequencies
        active.range_bin_size = ranges(2);
        active.time_axis = linspace(0,size(beat_frequncies,2)*active.pulse_duration,size(beat_frequncies,2));
      


%     % create CAF
%         range_doppler_caf = fftshift(fft(beat_frequncies,size(beat_frequncies,1)*2,2));
%         range_doppler_caf = transpose(20*log10(abs(range_doppler_caf./max(range_doppler_caf(:)))));
%         range_doppler_caf = range_doppler_caf(size(range_doppler_caf,2)/2:end,:);
%         figure
%         fig = imagesc(active.range_bins,active.time_axis,range_doppler_caf,[-120,0]); hold on;
%         ylabel('Time (s)');xlabel('Bistatic Range (m)');
%         c = colorbar;
%         c.Label.String='Norm Power (dB)'; 
% 



    % Derive range and time axis 
        active.n_range_bins = size(beat_frequncies,1);
        active.range_bins = 1:active.n_range_bins;
        active.fftfrequncies =fftfreq(size(beat_frequncies,1),1/(active.Fs/active.decimation_factor_actual)); % possible beat frequencies
        active.slope = active.Bw/active.pulse_duration;
        ranges = (active.fftfrequncies*C)/(2*active.slope); % calculate true range bin size    
        active.range_axis = ranges(1:(size(ranges,2)/2)); % truncate to only +ve beat frequencies
        active.range_bin_size = ranges(2)
        active.time_axis = linspace(0,size(beat_frequncies,2)*active.pulse_duration,size(beat_frequncies,2));
      
     % Plot RTI
        RTI_plot= transpose(20*log10(abs(beat_frequncies./max(beat_frequncies(:)))));
        figure
        fig = imagesc(active.range_bins,active.time_axis,RTI_plot,[-80,0]); hold on;
            ylabel('Time (s)');xlabel('Bistatic Range (m)');
            c = colorbar;
            c.Label.String='Norm Power (dB)'; 
            xlim([0 500])
            fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";

      % MTI Filtering 
        % Single Delay Line Filter 
            MTI_Data = zeros(size(beat_frequncies));
            active.range_bins = size(MTI_Data,1);
            % IIR Filter
                [b, a] = butter(12, 0.03, 'high');
                  for i=1:active.range_bins
                        MTI_Data(i,:) = filtfilt(b,a,beat_frequncies(i,:));
                  end
            
            
      % Plot RTI
        RTI_plot= transpose(20*log10(abs(MTI_Data./max(MTI_Data(:)))));
        figure
        fig = imagesc(1:active.range_bins,active.time_axis,RTI_plot,[-80,0]); hold on;
            ylabel('Time (s)');xlabel('Bistatic Range (m)');
            c = colorbar;
            c.Label.String='Norm Power (dB)'; 
%             xlim([0 200])


    phase_34 = angle(beat_frequncies(25,:));
    figure; plot(phase_34)


     % Plot series of pulse
        fig = figure
        ranges_2_plot = floor(linspace(1,100,4));
        ranges_2_plot = [1 100 200]
        for i = 1:length(ranges_2_plot)
            plot(active.range_axis*2,RTI_plot(ranges_2_plot(i),:));
            hold on
        end
            title("Single Pulse - " + Experiment_ID);
            xlim([0 2000])
            grid on; grid minor;
            ylabel('Relative Power (dB)')
            xlabel('Range (m)')  

            fig = figure
            plot(active.range_axis*2,RTI_plot(1,:),'black-'); hold on;
            plot(active.range_axis*2,RTI_plot(101,:),'r-');
            plot(active.range_axis*2,RTI_plot(201,:),'green-');
            xlim([0 500])
            grid on; grid minor;
            ylabel('Relative Power (dB)')
            xlabel('Bistatic Range (m)') 
            xline(108,'black--')
            leg = legend('Pulse 1 (t=0)','Pulse 101 (t=0.05s)','Pulse 201 (t=0.10s)','Breakthrough');
%             xline(275,'black--'),

%             title(leg,'Pulse Number:')
            fig_name = exp_dir + "Pulse_1_100_200" + Experiment_ID;
            saveas(fig,fig_name) 

            range_doppler_caf = fftshift(fft(MTI_Data,size(MTI_Data,1)*2,2));
            range_doppler_caf = transpose(20*log10(abs(range_doppler_caf./max(range_doppler_caf(:)))));
            figure
            fig = imagesc(active.range_bins,active.time_axis,range_doppler_caf,[-120,0]); hold on;
            ylabel('Time (s)');xlabel('Bistatic Range (m)');
            c = colorbar;
            c.Label.String='Norm Power (dB)'; 


 % Spectrogram 
         % Parameters
            r_start = 1;
            r_stop = 50;
            l_fft = 256;
            pad_factor = 2;
            overlap_factor = 0.90;
     
         % Plot Spectrogram pre-MTI filtering
            integrated_data = sum(beat_frequncies(r_start:r_stop,:));
            [spect,active.doppler_axis] = spectrogram(integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,active.PRF,'centered','yaxis');
            spect= 20*log10(abs(spect./max(spect(:))));
            figure
            fig = imagesc(active.time_axis,-active.doppler_axis,spect,[-50 0]);   
                ylim([-400 -100])
                xlim([-inf 5])
                c = colorbar
                c.Label.String='Norm Power (dB)'
                xlabel('Time (Sec)')
                % ylabel('Radial Velocity (mph)')   
                ylabel('Doppler Frequency (Hz)')  

end
                
 


%% phase correction 
% breakthrough = 5;
% target = 43;
breakthrough = 29;
target = 48;
power_noise = mean(abs(sum(beat_frequncies(500:560,:))))
power_breakthrough = mean(abs(beat_frequncies(breakthrough,:)));
power_target = mean(abs(beat_frequncies(target,:)));

snr_breakthrough = 20*log10(power_breakthrough)-20*log10(power_noise)
snr_target  = 20*log10(power_target )-20*log10(power_noise)

B_breakthrough = unwrap(angle(beat_frequncies(breakthrough,:)));
B_target = unwrap(angle(beat_frequncies(target,:)));

figure
plot(M_breakthrough,'blue-');hold on;
plot(M_target,'color','#4DBEEE');
plot(B_breakthrough,'red-');
plot(B_target,'Color','#D95319');
xlabel('Time (s)'); legend('Monostatic Breakthrough','Monostatic Target','Bistatic Breakthrough','Bistatic Target')
grid on; grid minor;

figure
plot(M_target - M_breakthrough,'blue-');hold on;
plot(B_target-B_breakthrough,'red-');
xlabel('Time (s)'); legend('Monostatic Breakthrough','Monostatic Target')
grid on; grid minor;


end