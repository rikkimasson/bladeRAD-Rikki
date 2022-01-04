
clear all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals/',...
        '/Users/piersbeasley/Documents/repos/bladeRAD/generic_scripts/matlab') % path to generic functions

%% Parameters - Configurable by User
save_directory1 = "/media/piers/data_drive/BladeRF_Experiments/7_Oct/hybrid/";
% save_directory1 = "/Volumes/data_drive/BladeRF_Experiments/7_Oct/hybrid/";

%create array with experiment names
dinfo = dir(save_directory1);
names_cell = {dinfo.name};
names_cell = names_cell(3:end);

process_fmcw = 0;
process_passive = 1;

for i = 17
    % load .mat file containing experiment parameters
        mat_file_name = save_directory1 + i + "/Hybrid Experimental Configuration.mat";
        load(mat_file_name);
   
        exp_dir = save_directory1 + i + '/';

%% FMCW Processing and Print RTI
   if process_fmcw == 1
       try     

                % load refsig for deramping
                    refsig = load_refsig(FMCW_Bw_M,FMCW_Fc,pulse_duration);    
                % load Signal, Mix and Dermap Signal  
                    file_location = exp_dir + 'active_' + Experiment_ID;
                    zero_padding = 2; % 1 = none; 2 = 100%
                    [max_range_actual,processed_signal] = deramp_and_decimate(file_location,FMCW_max_range,refsig,capture_duration,FMCW_number_pulses,FMCW_Fs,slope,zero_padding);
                    save(exp_dir + 'deramped_signal','processed_signal','max_range_actual')
                % Plot RTI
                    Range_axis = linspace(0,max_range_actual,size(processed_signal,1));
                    Range_bin = 1:size(processed_signal,1);
                    time_axis = linspace(0,size(processed_signal,2)*pulse_duration,size(processed_signal,2));
                    RTI_plot= transpose(10*log10(abs(processed_signal./max(processed_signal(:)))));
                    figure
                    fig = imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
                        ylabel('Time (Sec)')
                        xlabel('Range (m)')
                        %title("FMCW RTI - " + Experiment_ID)
                        xlim([0 200])
                        c = colorbar
                        c.Label.String='Norm Power (dB)'
                        fig_name = exp_dir + "RTI -" + Experiment_ID + ".jpg";
                        saveas(fig,fig_name,'jpeg')
                        fig_name = exp_dir + "RTI -" + Experiment_ID;
                        saveas(fig,fig_name)    


                  % Spectrogram 
                    single_rbin = 5;
                    start_rbin = 2;  
                    stop_rbin = 10;
                    spect_window = sum(processed_signal(start_rbin:stop_rbin,:));
                    l_fft = 512;
                    pad_factor = 4;
                    overlap_factor = 0.99;
            %         [spect,f] = spectrogram(processed_signal(single_rbin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
                    [spect,f] = spectrogram(spect_window,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');

                    % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
                    v=dop2speed(f,C/FMCW_Fc)*2.237;
                    spect= 10*log10(abs(spect./max(spect(:))));
                    figure
                    fig = imagesc(time_axis,-f,spect,[-30 0]);   
                        ylim([-100 100])
                        c = colorbar
                        c.Label.String='Norm Power (dB)'
                        xlabel('Time (Sec)')
                        % ylabel('Radial Velocity (mph)')   
                        ylabel('Doppler Frequency (Hz)')  
                        fig_title = "FMCW Spectrogram - R Bin: " + " - " + Experiment_ID;
                        title(fig_title);
                        fig_name = exp_dir + "FMCW Spectrogram_" + Experiment_ID + ".jpg";
                        saveas(fig,fig_name,'jpeg')
                        saveas(fig,fig_name)


            %% MTI Filtering 
                % Single Delay Line Filter 
                MTI_Data = zeros(size(processed_signal));
                      for i=2:FMCW_number_pulses
                            MTI_Data(:,i) = processed_signal(:,i)-processed_signal(:,i-1);
                      end

                % Plot MTI RTI      
                MTI_RTI_plot= transpose(10*log10(abs(MTI_Data./max(MTI_Data(:)))));
                figure
                fig = imagesc(Range_axis,time_axis,MTI_RTI_plot,[-20,0]);
                    xlim([1 200])
                    %ylim([0 0.0005])
                    grid on            
                    c = colorbar
                    c.Label.String='Norm Power (dB)'          
                    ylabel('Time (Sec)')
                    xlabel('Range (m)')   
                    fig_title = "Monostatic Single Delay Line MTI  RTI - Test " + Experiment_ID;
                    title(fig_title);
                    fig_name = exp_dir + "/MTI_RTI_" + Experiment_ID + ".jpg";
                    saveas(fig,fig_name,'jpeg')
                    fig_name = exp_dir + "/MTI_RTI_" + Experiment_ID;
                    saveas(fig,fig_name)
                    plot_signal = toc     


                %Plot MTI Spectrogram 
                single_rbin = 1;
                start_rbin = 2;  
                stop_rbin = 10;
                spect_window = sum(MTI_Data(start_rbin:stop_rbin,:));
                l_fft = 512;
                pad_factor = 4;
                overlap_factor = 0.90;
                %[spect,f] = spectrogram(MTI_Data(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
                [spect,f] = spectrogram(spect_window,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
                % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
                    v=dop2speed(f,C/FMCW_Fc)*2.237;
                    spect= 10*log10(abs(spect./max(spect(:))));
                    figure
                    fig = imagesc(time_axis,-f,spect,[-30 0]);
                    ylim([-100 100])
                    c = colorbar
                    c.Label.String='Norm Power (dB)'
                    xlabel('Time (Sec)')
                    % ylabel('Radial Velocity (mph)')   
                    ylabel('Doppler Frequency (Hz)')  
                    fig_title = "Monostatic Single Delay Line MTI Spectrogram - Test " + Experiment_ID;
                    title(fig_title);
                    fig_name = exp_dir + "/MTI_Spectrogram_" + Experiment_ID + ".jpg";
                    saveas(fig,fig_name,'jpeg')    
                    fig_name = exp_dir + "/MTI_Spectrogram_" + Experiment_ID;
                    saveas(fig,fig_name)   

    catch 
        "Error In FMCW Processing"
       end
   end
   
   
%% Passive Processing   
   if process_passive == 1
    try 

            % load signal and split ref and sur
                file_location = exp_dir + 'passive_' + Experiment_ID;
                [ref_channel, sur_channel]  = load_passive_data(file_location);
            % Plot time domain signals
                 figure
                 fig = subplot(2,1,1);
                    plot(real(ref_channel(1:samples_per_pulse)));
                    title("Ref channel time series");
                 hold on
                 subplot(2,1,2)
                     plot(real(sur_channel(1:samples_per_pulse)))
                     title("Sur channel time series");    
                     fig_name = exp_dir + "Time Domain Signals_" + Experiment_ID + ".jpg";
                     saveas(fig,fig_name,'jpeg')
                     fig_name = exp_dir + "Time Domain Signals_" + Experiment_ID;
                     saveas(fig,fig_name)
           % Batch process data and cross correlate  
                 seg_s = 5000; % number of segments per second - analagos to PRF.
                 seg_percent = 50;  % percentage of segment used for cross coreclation of 
                                    % survallance and reference. Will affect SNR dramatically.
                 [dec_ref_channel, cc_matrix] = passive_batch_process(ref_channel,sur_channel,seg_s,seg_percent,passive_Fs,passive_max_range,exp_dir);
                 save(exp_dir + 'passive_matrix','cc_matrix','seg_s','seg_percent','dec_ref_channel')
            % RTI Plot
                RTI_plot= transpose(10*log10(abs(cc_matrix./max(cc_matrix(:)))));
                Range_bin = linspace(0,passive_max_range,size(cc_matrix,1));
                
                time_axis = linspace(0,capture_duration,size(cc_matrix,2));
                figure
                fig = imagesc(Range_bin,time_axis,RTI_plot,[-50,0]);
                    % xlim([1 20])
                    %ylim([0 0.0005])
                    grid on            
                    c = colorbar
                    c.Label.String='Norm Power (dB)'
                    ylabel('Time (Sec)')
                    xlabel('Range (m)')   
                    fig_title = "Passive RTI - " + Experiment_ID;
                    title(fig_title);
                    fig_name = exp_dir + "Passive RTI_" + Experiment_ID + ".jpg";
                    saveas(fig,fig_name,'jpeg')
                    fig_name = exp_dir + "Passive RTI_" + Experiment_ID;
                    saveas(fig,fig_name)

              % CAF of entire capture
                f_axis = linspace(-seg_s/2,seg_s/2,size(cc_matrix,2));
                t_cc_matrix = transpose(cc_matrix);
                CAF = fftshift(fft(t_cc_matrix,size(t_cc_matrix,1),1),1);
                figure
                imagesc(Range_bin,f_axis,10*log10(abs(CAF./max(CAF(:)))),[-50 1]); 
                    ylim([-500 500])     
                    % xlim([1 20])
                    colorbar
                    ylabel('Doppler Shift (Hz)')
                    xlabel('Range Bin')  
                    title("CAF for entire capture" + Experiment_ID)


             % Spectrogram 
                single_rbin = 5;
                start_rbin = 2;  
                stop_rbin = 10;
                spect_window = sum(cc_matrix(start_rbin:stop_rbin,:));
                l_fft = 512;
                pad_factor = 4;
                overlap_factor = 0.90;
        %         [spect,f] = spectrogram(single_rbin,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,seg_s,'centered','yaxis');
                [spect,f] = spectrogram(spect_window,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,seg_s,'centered','yaxis');
                % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
                v=dop2speed(f,C/passive_Fc)*2.237;
                spect= 10*log10(abs(spect./max(spect(:))));
                figure
                fig = imagesc(time_axis,f,spect,[-30 0]);   
                    ylim([-200 200])
                    c = colorbar
                    c.Label.String='Norm Power (dB)'
                    xlabel('Time (Sec)')
                    % ylabel('Radial Velocity (mph)')   
                    ylabel('Doppler Frequency (Hz)')  
                    fig_title = "Passive Spectrogram - R Bin:" + Experiment_ID;
                    title(fig_title);
                    fig_name = exp_dir + "Passive Spectrogram_" + Experiment_ID + ".jpg";
                    saveas(fig,fig_name,'jpeg')
                    fig_name = exp_dir + "Passive Spectrogram_" + Experiment_ID;
                    saveas(fig,fig_name)
    catch 
        "Error In Passive Processing"
    end
   end
   
   
end

%% Save Passive Files as int array 



inter_sur = interleave(real(sur_channel),imag(sur_channel));
inter_sur = inter_sur * 1e4;


inter_ref = interleave(real(ref_channel),imag(ref_channel));
inter_ref = inter_ref * 1e4;
