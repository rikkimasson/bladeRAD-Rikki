clear all
addpath('/home/piers/repos/bladeRAD/generic_scripts/matlab',...
        '/home/piers/repos/bladeRAD/generic_scripts',...
        '/home/piers/repos/bladeRAD/generic_scripts/ref_signals/') % path to generic functions

%% Parameters - Configurable by User

% Capture parameters 
Experiment_ID = 0005;       % Expeiment Name
capture_duration = 10;    % capture duration
Bw = 15e6;               % Sample Rate of SDR per I & Q (in reality Fs is double this)
save_directory = "/home/piers/Documents/Captures/"; % each experiment will save as a new folder in this directory
passive_max_range = 1000; %max range to cross-correlate to

% Radar Parameters 
Fc = 530e6;   % Central RF    
Ref_gain = 47;
Sur_gain = 60;
Pass_SDR = 3;   % SDR to use for Passive Radar - labelled on RFIC Cover and bladeRAD Facia Panel

% Parameters not configurable by user 
    C = physconst('LightSpeed');
    Fs = Bw;
    sample_duration = 1/Fs;
    number_cap_samps = 2*(capture_duration/sample_duration)
    RF_freq = Fc/1e6;   % RF in MHz 
    Bw_M = Bw/1e6;      % BW in MHz
    file_size_MBytes = (number_cap_samps * 16)*2/(8*1e6) 

    
%% Setup Radar
    % 1 'set clock_sel external'; 2 'set clock_out enable; 3 'set clock_ref enable'

    % Setup Passive SDR 
   [trig_flag, passive_command] = create_shell_command(Experiment_ID,...
                                   number_cap_samps,... 
                                   0,...
                                   0,...
                                   0,...
                                   Ref_gain,...
                                   Sur_gain,...
                                   RF_freq,...
                                   Bw_M,...
                                   Pass_SDR,...
                                   'master',...
                                   3,...
                                   'pass');
                               
    %passive_command = tx_command + "&"; % uncomment for non-blocking system command execution                    
    status = system(passive_command);

%% Save Raw Data and create header file for directory 
    exp_dir = save_directory + Experiment_ID + '/';
    make_dir = 'mkdir ' + exp_dir;
    system(make_dir); % Blocking system command execution
    move_file = 'mv /tmp/passive_' + string(Experiment_ID) + '.sc16q11 ' + exp_dir;
    rtn = system(move_file);
    if rtn == 0
        "Rx Data Copyied to Save directory"
    else 
        "Rx Copy Failed"
        return
    end
    save(exp_dir + 'Passive Experimental Configuration') 


 %% Passive Processing
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
         seg_s = 1000; % number of segments per second - analagos to PRF.
         seg_percent = 90;  % percentage of segment used for cross coreclation of 
                            % survallance and reference. Will affect SNR dramatically.
         cc_matrix = passive_batch_process(ref_channel,sur_channel,seg_s,seg_percent,Fs,passive_max_range,exp_dir);
         save(exp_dir + 'passive_matrix','cc_matrix')
    % RTI Plot
        RTI_plot= transpose(10*log10(abs(cc_matrix./max(cc_matrix(:)))));
        Range_bin = linspace(0,passive_max_range,size(cc_matrix,1));
        time_axis = linspace(0,capture_duration,size(cc_matrix,2));
        figure
        fig = imagesc(Range_bin,time_axis,RTI_plot,[-50,0]);
            % xlim([1 20])
            %ylim([0 0.0005])
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
        f_axis = linspace(-seg_s/2,seg_s/2,size(cc_matrix,2));
        t_cc_matrix = transpose(cc_matrix);
        CAF = fftshift(fft(t_cc_matrix,size(t_cc_matrix,1),1),1);
        figure
        imagesc(Range_bin,f_axis,10*log10(abs(CAF./max(CAF(:)))),[-100 1]); 
            ylim([-500 500])     
            % xlim([1 20])
            colorbar
            ylabel('Doppler Shift (Hz)')
            xlabel('Range Bin')  
            title("CAF for entire capture" + Experiment_ID)
            fig_name = exp_dir + "CAF for entire capture_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)
        

     % Spectrogram 
        r_bin = 1;
        l_fft = 2056;
        pad_factor = 1;
        overlap_factor = 0.99;
        [spect,f] = spectrogram(cc_matrix(r_bin,:),l_fft,round(l_fft*overlap_factor),l_fft*pad_factor,seg_s,'centered','yaxis');
        % spect(pad_factor*l_fft/2-1:pad_factor*l_fft/2+1,:) = 0;
        v=dop2speed(f,C/Fc)*2.237;
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
        fig = imagesc(time_axis,f,spect,[-30 0]);   
            ylim([-600 600])
            colorbar
            xlabel('Time (Sec)')
            % ylabel('Radial Velocity (mph)')   
            ylabel('Doppler Frequency (Hz)')  
            fig_title = "Passive Spectrogram - R Bin: " + r_bin + " - " + Experiment_ID;
            title(fig_title);
            fig_name = exp_dir + "Passive Spectrogram_" + Experiment_ID + ".jpg";
            saveas(fig,fig_name,'jpeg')
            saveas(fig,fig_name)



         