clear all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 

% save_directory1 = "/home/piers/Documents/Captures/Pulse_Doppler/"; % rach experiment will save as a new folder in this directory
save_directory1 = "/media/piers/T7/15_09_2022_farm/pulse_doppler/"; % each experiment will save as a new folder in this directory

% create array with experiment names
dinfo = dir(save_directory1);
names_cell = {dinfo.name};
names_cell = names_cell(3:end);

for i = 65:66
 Experiment_ID = i;

 %% load .mat file containing experiment parameters
    mat_file_name = save_directory1 + i + "/Experimental Configuration.mat";
    load(mat_file_name);  

 %% Load Reference Deramp Signal
    exp_dir = save_directory1 + Experiment_ID + '/';
    refsig = load_refsig(Bw_M,Fs,Fc,pulse_duration);
    
%% Load Signals
    file_location = exp_dir + 'active_' + Experiment_ID + '.sc16q11';
    raw_data = load_sc16q11(file_location);    
    
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
zero_padding = 1;

% Match Filter 
matched_filter = refsig;
matched_filter_fft = conj(fft(matched_filter,(samples_per_PRI*zero_padding))); 
processed_signal = zeros(size(pulse_matrix,2),size(matched_filter_fft,1));
for i = 1:number_pulses
    pri = fft(pulse_matrix(:,i),size(matched_filter_fft,1));
    processed_signal(i,:) = ifft(pri.*matched_filter_fft); 
end

%% PLot RTI
    Range_bins = 0:size(processed_signal,2)-1;
    Range_bin_size = 1/(Fs*2) *C;
    Range_axis = Range_bins * Range_bin_size;
    time_axis = linspace(0,capture_duration,size(processed_signal,1));
    RTI_plot=10*log10(abs(processed_signal./max(max(processed_signal))));
    figure
    fig = imagesc(Range_axis,time_axis,RTI_plot,[-50,0]); 
        xlim([0 100]);
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
        plot(Range_axis,RTI_plot(i,:));
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
                    MTI_Data(i,1:200) = processed_signal(i,1:200)-processed_signal(i-1,1:200);
              end        
        


  %% PLot MTI RTI
    RTI_plot=10*log10(abs(MTI_Data./max(max(MTI_Data))));
    figure
    fig = imagesc(Range_axis,time_axis,RTI_plot,[-50,0]); 
        xlim([0 1000]);
        grid on            
        colorbar
        ylabel('Time (Sec)')
        xlabel('Range')      
        title("Pulse Doppler - Monostatic RTI - " + Experiment_ID);
        fig_name = exp_dir + "Monostatic_RTI_" + Experiment_ID + ".jpg";
        saveas(fig,fig_name,'jpeg')
    
%% Spectrogram rom
        r_start = 5;
        r_stop = 15;
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



end