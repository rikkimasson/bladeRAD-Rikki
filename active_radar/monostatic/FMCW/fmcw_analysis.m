% clear all
close all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
% save_directory1 = "/media/piers/data_drive/BladeRF_Experiments/7_Oct/FMCW/"; % each experiment will save as a new folder in this directory
save_directory1 = "/media/piers/T7/experiments lab 28_06_22/multistatic experiments/Lab Multistatic Experiments N1/";

% create array with experiment names
dinfo = dir(save_directory1);
names_cell = {dinfo.name};
names_cell = names_cell(3:end);
i=37
for i = 37
%     
%     % load .mat file containing experiment parameters
%     mat_file_name = save_directory1 + i + "/Experimental Configuration.mat";
%     load(mat_file_name);
   
    %% Load Reference Deramp Signal
        refsig = load_refsig(Bw_M,Fs,Fc,pulse_duration);

    %% Load Signal, Mix and Dermap Signal 
        exp_dir = save_directory1 + i + '/';
        file_location = exp_dir + 'active_' + Experiment_ID;
        [max_range_actual,decimated_signal] = deramp_and_decimate(file_location,max_range,refsig,capture_duration,number_pulses,Fs,slope,true);
        save(exp_dir + 'deramped_signal','decimated_signal','max_range_actual')

    %% Window and FFT Signal 
    % window signal
        w = window('hann',size(decimated_signal,1));
        windowed_signal = decimated_signal.*w;
    % fft signal
        zero_padding = 1;
        processed_signal = fft(windowed_signal,size(windowed_signal,1)*zero_padding);
    
    
    %% MTI Filtering 
        % Single Delay Line Filter 
        MTI_Data = zeros(size(processed_signal));
              for i=2:number_pulses
                    MTI_Data(:,i) = processed_signal(:,i)-processed_signal(:,i-1);
              end
        
     %% Plot RTI

        Range_axis = linspace(0,max_range_actual,size(processed_signal,1));
        Range_bin = 1:size(processed_signal,1);
        time_axis = linspace(0,size(processed_signal,2)*pulse_duration,size(processed_signal,2));
        RTI_plot= transpose(10*log10(abs(processed_signal./max(processed_signal(:)))));
        figure
        fig = imagesc(Range_axis,time_axis,RTI_plot,[-50,0]);   
            xlim([-inf 100])
            grid on            
            colorbar
            ylabel('Time (Sec)')
            xlabel('Range (m)')   
%             title("FMCW RTI - " + Experiment_ID);
            c = colorbar
            c.Label.String='Normalised Power (dB)';
            fig_name = exp_dir + "FMCW_RTI_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg') 



        figure
        fig = plot(Range_axis,RTI_plot(10,:));
            title("Single Pulse - " + Experiment_ID);
            xlim([0 100])
            grid on
            ylabel('Relative Power (dB)')
            xlabel('Range (m)')  
            fig_name = exp_dir + "Single_Pulse_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg') 
            fig_name = exp_dir + "Single_Pulse_" + Experiment_ID;
            saveas(fig,fig_name) 

    % Plot Spectrogram 
        single_rbin = 1;
        start_rbin = 2;  
        stop_rbin = 10;
        spect_window = sum(processed_signal(start_rbin:stop_rbin,:));
        l_fft = 512;
        pad_factor = 1;
        overlap_factor = 0.99;
        %[spect,f] = spectrogram(processed_signal(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
        [spect,f] = spectrogram(spect_window,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
            v=dop2speed(f,C/Fc)*2.237;
            spect= 10*log10(abs(spect./max(spect(:))));
            figure
            fig = imagesc(time_axis,f,spect,[-50 0]);
            ylim([-100 100])
            colorbar
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "Monostatic Spectrogram - Test " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "/Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')        
            fig_name = exp_dir + "/Spectrogram_" + Experiment_ID;
            saveas(fig,fig_name) 

    %% MTI Filtering 
        % Single Delay Line Filter 
        MTI_Data = zeros(size(processed_signal));
              for i=2:number_pulses
                    MTI_Data(:,i) = processed_signal(:,i)-processed_signal(:,i-1);
              end

        % Plot MTI RTI      
        MTI_RTI_plot= transpose(10*log10(abs(MTI_Data./max(MTI_Data(:)))));
        figure
        fig = imagesc(Range_axis,time_axis,MTI_RTI_plot,[-30,0]);
            xlim([1 100])
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
        spect_window = sum(MTI_Data(start_rbin:stop_rbin,:));
        overlap_factor = 0.99;
        %[spect,f] = spectrogram(MTI_Data(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
        [spect,f] = spectrogram(spect_window,l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,PRF,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
            v=dop2speed(f,C/Fc)*2.237;
            spect= 10*log10(abs(spect./max(spect(:))));
            figure
            fig = imagesc(time_axis,f,spect,[-30 0]);
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


    %     %% Coherent integration 
    %     compressed_data = sum(Dec_Deramped,2);
    %     figure
    %     plot(Range_axis,abs(compressed_data))
    %     xlim([0 1000])

end

%% Phase Analysis
node1_phase = angle(processed_signal(4,:));

conv_factor = 1/(2*pi*Fc)
y = node0_phase * conv_factor;
x = time_axis;
coefficients = polyfit(x,y,1)
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);

figure
plot(time_axis,node0_phase*conv_factor+3e-11)
hold on
plot(xFit,yFit+3e-11,'Color','black')
hold on 
plot(time_axis,node1_phase*conv_factor,'Color','#D95319')
ylabel('Phase Time (s)')
grid on 
grid minor
xlabel('Time (s)')
legend('Monostatic Phase','Regression Line','Bistatic Phase')

nfft = 2^15
noverlap = 0.8;
[pxx,f] = pwelch(node0_phase,hann(nfft),0,nfft,'psd',PRF);
[pxx1,f1] = pwelch(node1_phase,hann(nfft),0,nfft,'psd',PRF);


figure
semilogx(f,10*log10(pxx./max(pxx)))
hold on 
semilogx(f1,10*log10(pxx1./max(pxx1)))
grid on
grid minor
xlabel('Frequency (Hz)')
ylabel('Phase Noise (dBc/Hz)')
legend('Monostatic Phase','Bistatic Phase')
xlim([-inf 2500])



