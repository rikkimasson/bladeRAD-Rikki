% clear all
% % close all
addpath('~/repos/bladeRAD/generic_scripts/matlab',...
        '~/repos/bladeRAD/generic_scripts',...
        '~/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

process_active_a = true;
process_passive_a = true;
 
%% create array with experiment names
dinfo = dir(save_directory1);
names_cell = {dinfo.name};
names_cell = names_cell(3:end);

for i=2
    
%% load .mat file containing experiment parameters
        mat_file_name = save_directory1 + i + "/Experimental Configuration.mat";
        load(mat_file_name);
   
        exp_dir = save_directory1 + i + '/';

%% FMCW Processing and Print RTI
if process_active_a == true
    % load refsig for deramping
        refsig = load_refsig(active.Bw_M,active.Fs,active.Fc,active.pulse_duration);    
    % load Signal, Mix and Dermap Signal  
        file_location = exp_dir + 'active_' + Experiment_ID;
        lp_filter = true;
        [max_range_actual,deramped_signal] = deramp_and_decimate(file_location,active.max_range,refsig,capture_duration,active.number_pulses,active.Fs,active.slope,lp_filter);
        save(exp_dir + 'deramped_signal','deramped_signal')
    
    % Window and FFT Signal 
    % window signal
        w = window('hamming',size(deramped_signal,1));
        windowed_signal = deramped_signal.*w;
    % fft signal
        zero_padding = 1; % 1 = none; 2 = 100%
        processed_signal = fft(windowed_signal,size(windowed_signal,1)*zero_padding);

    % MTI Filtering 
        % Single Delay Line Filter 
        MTI_Data = zeros(size(processed_signal));
              for i=2:active.number_pulses
                    MTI_Data(:,i) = processed_signal(:,i)-processed_signal(:,i-1);
              end
        
      % Plot RTI
        Range_axis = linspace(0,max_range_actual,size(processed_signal,1));
        Range_bin = 1:size(processed_signal,1);
        time_axis = linspace(0,size(processed_signal,2)*active.pulse_duration,size(processed_signal,2));
        RTI_plot= transpose(10*log10(abs(processed_signal./max(processed_signal(:)))));
        figure
        fig = imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
            ylabel('Time (Sec)')
            xlabel('Range (m)')
            title("FMCW RTI - " + Experiment_ID)
            xlim([0 200])
            c = colorbar
            c.Label.String='Norm Power (dB)'
            fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)    
    
         MTI_RTI_plot= transpose(10*log10(abs(MTI_Data./max(MTI_Data(:)))));
                figure
                fig = imagesc(Range_axis,time_axis,MTI_RTI_plot,[-50,0]);
                    xlim([1 200])
                    %ylim([0 0.0005])
                    grid on            
                    c = colorbar
                    c.Label.String='Norm Power (dB)'          
                    ylabel('Time (Sec)')
                    xlabel('Range (m)')   
                    fig_title = "MTI  RTI - " + Experiment_ID;
                    title(fig_title);
                    fig_name = exp_dir + "/MTI_RTI_" + Experiment_ID + ".jpg";
                    saveas(fig,fig_name,'jpeg')
                    fig_name = exp_dir + "/MTI_RTI_" + Experiment_ID;
                    saveas(fig,fig_name)


    fig = figure
    ranges_2_plot = floor(linspace(1,active.number_pulses,5));
    for i = ranges_2_plot
        plot(RTI_plot(i,:));
        hold on
    end
        title("Single Pulse - " + Experiment_ID);
        xlim([0 200])
        grid on
        ylabel('Relative Power (dB)')
        xlabel('Range (m)')  
        fig_name = exp_dir + "Single_Pulse" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg') 

            
     
      % Spectrogram 
        r_start = 2;
        r_stop = 10;
        l_fft = 256;
        pad_factor = 1;
        overlap_factor = 0.90;
        
        integrated_data = sum(processed_signal(r_start:r_stop,:));
        [spect,f] = spectrogram(integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,active.PRF,'centered','yaxis');
        v=dop2speed(f,C/passive.Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,-f,spect,[-50 0]);   
            ylim([-500 500])
            c = colorbar
            c.Label.String='Norm Power (dB)'
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "FMCW Spectrogram -" + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID;
            saveas(fig,fig_name)
            
            MTI_integrated_data = sum(MTI_Data(r_start:r_stop,:));
            [spect,f] = spectrogram(MTI_integrated_data,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,active.PRF,'centered','yaxis');
            v=dop2speed(f,C/passive.Fc)*2.237;
            spect= 10*log10(abs(spect./max(spect(:))));
            figure
            fig = imagesc(time_axis,-f,spect,[-50 0]);   
                ylim([-500 500])
                c = colorbar
                c.Label.String='Norm Power (dB)'
                xlabel('Time (Sec)')
                % ylabel('Radial Velocity (mph)')   
                ylabel('Doppler Frequency (Hz)')  
                fig_title = "MTI FMCW Spectrogram - " + Experiment_ID;
                title(fig_title);
                fig_name = exp_dir + "MTI FMCW Spectrogram_" + Experiment_ID + ".jpg";
                saveas(fig,fig_name,'jpeg')
                fig_name = exp_dir + "MTI FMCW Spectrogram_" + Experiment_ID;
                saveas(fig,fig_name)




%% Process Active data into Range-Doppler Slices
           active.cpi = 0.5; % cohernet proccessing interval (s)
           active.cpi_overlap = 0.5; % overlap between CPIs (watch this - too large will cause slow exceution)
           active.zero_padding = 1;
           active.doppler_window = 'blackmanharris';
           active.dynamic_range = 50;
           active.velocity_conv = C*(((1/C)/(active.Fc/C)));
                 active.range_bins = size(MTI_Data,1);
                 active.max_range_actual = max_range_actual; %1000;
  
          [active.number_cpi,....
           active.pulses_per_cpi,...
           active.range_doppler_slices] = rangeDopplerSlice(MTI_Data,active.cpi,active.PRF,...
                                                                 active.cpi_overlap,...
                                                                 active.zero_padding,...
                                                                 active.doppler_window);   
           
           % Determine Active Data Range and Doppler Axis                                              
             active.doppler_bins = active.pulses_per_cpi*active.zero_padding;
             active.doppler_axis = linspace(-active.PRF/2,active.PRF/2,active.doppler_bins);
             active.doppler_velocity_axis = active.doppler_axis * active.velocity_conv;
             active.range_axis = linspace(0,active.max_range_actual,active.range_bins);
             active.range_bins_axis = 1:active.range_bins;

           % Create video of range-Doppler slices
             video_name = exp_dir + "active_range-Doppler_0ZP_hann_wind_Dopppler_Exp_" + Experiment_ID + ".avi";
             %video_name = "range-Doppler_log_Exp_" + Experiment_ID + ".avi";       
             video_title = "Active Radar Capture";
             dynamic_range = +50;
             max_range = 100;
             max_doppler = 200;
             frame_rate = 1/(capture_duration/active.number_cpi);    

             createVideo(active.range_doppler_slices,frame_rate,...
                         active.range_axis,max_range,...
                         -active.doppler_axis,max_doppler,...
                         dynamic_range,video_name,video_title);




                
end

 %% Passive Processing
 if process_passive_a == true
    passive.max_range = 70
    passive.seg_s = 5000; % number of segments per second - analagos to PRF.
    passive.seg_percent = 100;  % percentage of segment used for cross coreclation of 
                            % survallance and reference. Will affect SNR dramatically.
    % load signal and split ref and sur
        file_location = exp_dir + 'passive_' + Experiment_ID;
        [ref_channel, sur_channel]  = load_passive_data(file_location);
    % Plot time domain signals
         figure
         fig = subplot(2,1,1);
            plot(real(ref_channel(1:4000000)));
            title("Ref channel time series");
         hold on
         subplot(2,1,2)
             plot(real(sur_channel(1:4000000)))
             title("Sur channel time series");    
             fig_name = exp_dir + "Time Domain Signals_" + Experiment_ID + ".jpg";
             saveas(fig,fig_name,'jpeg')
   % Batch process data and cross correlate  
         [ref_matrix ,self_ambg_matrix, cc_matrix] = passive_batch_process(ref_channel,sur_channel,passive.seg_s,passive.seg_percent,passive.Fs,passive.max_range,exp_dir);
         save(exp_dir + 'passive_matrix','cc_matrix')
    % RTI Plot
        RTI_plot= transpose(10*log10(abs(cc_matrix./max(cc_matrix(:)))));
        Range_bin = linspace(0,passive.max_range,size(cc_matrix,1));
        time_axis = linspace(0,capture_duration,size(cc_matrix,2));
        figure
        fig = imagesc(Range_bin,time_axis,RTI_plot,[-50,0]);
            xlim([1 20])
%             ylim([0 0.0005])
            grid on            
            colorbar
            ylabel('Time (Sec)')
            xlabel('Range Bin')   
            fig_title = "Passive RTI - " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "Passive RTI_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)

      % CAF of entire capture
        f_axis = linspace(-passive.seg_s/2,passive.seg_s/2,size(cc_matrix,2));
        t_cc_matrix = transpose(cc_matrix);
        CAF = fftshift(fft(t_cc_matrix,size(t_cc_matrix,1),1),1);
        figure
        imagesc(Range_bin,f_axis,10*log10(abs(CAF./max(CAF(:)))),[-50 1]); 
            ylim([-500 500])     
            xlim([1 20])
            c = colorbar
            c.Label.String='Norm Power (dB)'
            ylabel('Doppler Shift (Hz)')
            xlabel('Range Bin')  
            title("CAF for entire capture" + Experiment_ID)
            fig_name = exp_dir + "CAF for entire capture_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)
        

     % Spectrogram 
        int_bins = sum(cc_matrix(2:10,:),1);
        r_bin = 1;
        l_fft = 1024;
        pad_factor = 1;
        overlap_factor = 0.99;
        [spect,f] = spectrogram(int_bins,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,passive.seg_s,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/passive.Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-30 0]);   
            ylim([-500 +500])
            c = colorbar
            c.Label.String='Norm Power (dB)'
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "Passive Spectrogram :- " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "Passive Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)


%% Proccess Passive data into Range-Doppler Slices
           passive.PRF = passive.seg_s; %seg_s
           passive.cpi = 0.5; % cohernet proccessing interval (s)
           passive.cpi_overlap = 0.5; % overlap between CPIs (watch this - too large will cause slow exceution)
           passive.doppler_window = 'hann';
           passive.zero_padding = 1;
           passive.dynamic_range = +50;
           passive.max_range_m = passive.max_range*(1/passive.Fs)/passive.zero_padding * C;      
     
                [passive.number_cpi,... 
                 passive.pulses_per_cpi,...
                 passive.range_doppler_slices] = rangeDopplerSlice(cc_matrix,passive.cpi,passive.PRF,...
                                                                   passive.cpi_overlap,...
                                                                   passive.zero_padding,...
                                                                   passive.doppler_window);  
                
                 [~,~,passive.self_ambg_slices] = rangeDopplerSlice(self_ambg_matrix,passive.cpi,passive.PRF,...
                                                                   passive.cpi_overlap,...
                                                                   passive.zero_padding,...
                                                                   passive.doppler_window);  
            
            passive.cpi_stride = round(passive.pulses_per_cpi*(1-passive.cpi_overlap)); % number of pulses to stride each for next CPI
            passive.velocity_conv = C*(((1/C)/(passive.Fc/C)));
            passive.range_bins = size(cc_matrix,1);
            passive.doppler_bins = passive.pulses_per_cpi*passive.zero_padding+1;
            passive.doppler_axis = linspace(-passive.PRF/2,passive.PRF/2,passive.doppler_bins);
            passive.doppler_velocity_axis = passive.doppler_axis*passive.velocity_conv;
            passive.range_axis = linspace(0,passive.max_range_m,passive.range_bins);


            % create video of  range-Doppler slices
             video_name = exp_dir + "Range_Doppler_Slices" + Experiment_ID + ".avi";
             %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";       
             video_title = "Passive Pre-DSI";
             dynamic_range = +50;
             max_range = 200;
             max_doppler = 100;
             frame_rate = 1/(capture_duration/passive.number_cpi);    
             createVideo(passive.range_doppler_slices,frame_rate,...
                         passive.range_axis,max_range,...
                         passive.doppler_axis,max_doppler,...
                         dynamic_range,video_name,video_title);
        
                 
%% Direct Signal Interference Cancellation
   % set DSI cancellation parameters
      p = 0.999;        % subtraction parameter - P must be a positive 
                        % integer less than one to avoid unwanted discontinuities
                        % arising from zero values in the range–Doppler surface.
     threshold = 0.005; % cutoff threshold parameter
     max_iterations = 100; % maximum number of itterations DSI is CLEANed from CAF slice
     number_rbins = size(passive.range_doppler_slices{1},2);
  
   % perform CLEAN based DSI Cancellation   
     passive.CLEANed_range_doppler_slices = CLEAN_DSI(passive.range_doppler_slices,...
                                               passive.self_ambg_slices,...
                                               ref_matrix,...
                                               passive.number_cpi,...
                                               passive.pulses_per_cpi,...
                                               passive.cpi_stride,...
                                               passive.zero_padding,...
                                               number_rbins,...
                                               max_iterations,...
                                               threshold,p,...
                                               passive.range_axis,passive.doppler_axis);

close all
  % create video of CLEANed range-Doppler slices
     video_name = exp_dir + "CLEANed_range-Doppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";
     %video_name = "passive_RangeDoppler_CLEANed_log_Exp_" + Experiment_ID + ".avi";       
     video_title = "CLEANing Passive Radar Capture";
     dynamic_range = 50;
     max_range = 200;
     max_doppler = 100;
     frame_rate = 1/(capture_duration/passive.number_cpi);    
     createVideo(passive.CLEANed_range_doppler_slices,frame_rate,...
                 passive.range_axis,max_range,...
                 passive.doppler_axis,max_doppler,...
                 dynamic_range,video_name,video_title);

 end

end
