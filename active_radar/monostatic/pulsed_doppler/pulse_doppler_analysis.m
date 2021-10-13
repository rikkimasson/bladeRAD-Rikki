clear all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 

% save_directory1 = "/home/piers/Documents/Captures/Pulse_Doppler/"; % rach experiment will save as a new folder in this directory
save_directory1 = "/media/piers/data_drive/BladeRF_Experiments/7_Oct/pulse_doppler/"; % each experiment will save as a new folder in this directory

% create array with experiment names
dinfo = dir(save_directory1);
names_cell = {dinfo.name};
names_cell = names_cell(3:end);

for i = 4
    % load .mat file containing experiment parameters
        mat_file_name = save_directory1 + i + "/FMCW Experimental Configuration.mat";
        load(mat_file_name);
        
        exp_dir = save_directory1 + i + '/';

          %% Load Reference Deramp Signal
            refsig = load_refsig(Bw_M,Fc,pulse_duration);

        %% Load Signals
            file_location = exp_dir + 'active_' + Experiment_ID + '.sc16q11';
            raw_data = load_sc16q11(file_location);


        %% Reshape array into matrix of pulses
        pulse_matrix = reshape(raw_data,[length(raw_data)/number_pulses,number_pulses]); %reshape array to individual pulses
        %figure
        %spectrogram(pulse_matrix(:,3),128,100,100,Fs,'centered','yaxis')
        figure
        fig = plot(abs(real(pulse_matrix(:,3))));
            ylabel('ADC Value (0-1)')
            xlabel('Samples')      
            title("Pulse Doppler - Pulse Time Series - " + Experiment_ID);
            fig_name = exp_dir + "Pulse_Time_Series_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            fig_name = exp_dir + "Pulse_Time_Series_" + Experiment_ID;
            saveas(fig,fig_name)


        %% Range Limmited Signal Processing 
        zero_padding = 1;
        % Match Filter 
        matched_filter = refsig;
        matched_filter_fft = conj(fft(matched_filter,(samples_per_PRI*zero_padding))); 
        radar_matrix = zeros(size(pulse_matrix,2),size(matched_filter_fft,1));

        for i = 1:number_pulses
            pri = fft(pulse_matrix(:,i),size(matched_filter_fft,1));
            radar_matrix(i,:) = ifft(pri.*matched_filter_fft); 
        end
            radar_matrix = radar_matrix(:,1:2000);
            Range_bin = 1:size(radar_matrix,2);
            Range = linspace(0,R_Max*2,size(radar_matrix,2));
            time_axis = linspace(0,capture_duration,size(radar_matrix,1));
            RTI_plot=10*log10(abs(radar_matrix./max(radar_matrix(:))));
            figure
            fig = imagesc(Range,time_axis,RTI_plot,[-50,0]); 
%                 xlim([0 100]);
                grid on            
                colorbar
                ylabel('Time (Sec)')
                xlabel('Range')      
                title("Pulse Doppler - Monostatic RTI - " + Experiment_ID);
                fig_name = exp_dir + "Monostatic_RTI_" + Experiment_ID + ".jpg";
                saveas(fig,fig_name,'jpeg')
                fig_name = exp_dir + "Monostatic_RTI_" + Experiment_ID;
                saveas(fig,fig_name)                
            figure
            plot(RTI_plot(100,:));
                title("Single Pulse - " + Experiment_ID);
                xlim([0 100])
                grid on
                ylabel('Relative Power (dB)')
                xlabel('Range (m)')  
                fig_name = exp_dir + "Single_Pulse" + Experiment_ID + ".jpg";
                saveas(fig,fig_name,'jpeg') 
                fig_name = exp_dir + "Single_Pulse" + Experiment_ID;
                saveas(fig,fig_name) 
                
                %% MTI Filtering 
                % Single Delay Line Filter 
                MTI_Data = zeros(size(radar_matrix));
                      for i=2:number_pulses
                            MTI_Data(i,:) = radar_matrix(i,:)-radar_matrix(i-1,:);
                      end

                % Plot MTI RTI      
                MTI_RTI_plot= 10*log10(abs(MTI_Data./max(MTI_Data(:))));
                figure
                fig = imagesc(Range_bin,time_axis,MTI_RTI_plot,[-20,0]);
                    xlim([1 50])
                    %ylim([0 0.0005])
                    grid on            
                    colorbar
                    ylabel('Time (Sec)')
                    xlabel('Range Bin')   
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
                MTI_Data_t = transpose(MTI_Data);
                spect_window = sum(MTI_Data_t(start_rbin:stop_rbin,:),1);
                l_fft = 512;
                pad_factor = 4;
                overlap_factor = 0.90;
                %[spect,f] = spectrogram(MTI_Data(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
                [spect,f] = spectrogram(spect_window,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
                % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
                    spect= 10*log10(abs(spect./max(spect(:))));
                    figure
                    fig = imagesc(time_axis,f,spect,[-25 0]);
                    ylim([-100 100])
                    colorbar
                    xlabel('Time (Sec)')
                    % ylabel('Radial Velocity (mph)')   
                    ylabel('Doppler Frequency (Hz)')  
                    fig_title = "Monostatic Single Delay Line MTI Spectrogram - Test " + Experiment_ID;
                    title(fig_title);
                    fig_name = exp_dir + "/MTI_Spectrogram_" + Experiment_ID + ".jpg";
                    saveas(fig,fig_name,'jpeg')    
                    fig_name = exp_dir + "/MTI_Spectrogram_" + Experiment_ID;
                    saveas(fig,fig_name)   

                
                

end