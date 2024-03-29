function [actual_max_range,decimated_signal,decimation_factor_actual] = deramp_and_decimate(file_location,max_range,refsig,capture_duration,number_pulses,Fs,slope,lp_filter)
% MIX_AND_DERAMP - for bladeRAD FMCW Signal Processing  

%   1. Load receive file - location passed in as 'file_location'
%   2. Reshape receive channel into radar matrix of PRI
%   3. Deramp signal using a pre-recorded chirp passed in as'refsig'.
%   4. Decimate signal - low pass and decimate to range passed in as 'max_range'


C = physconst('LightSpeed');

if capture_duration < 61
     
    %% load rx file
            rx_file= file_location + ".sc16q11";
            raw_data = load_sc16q11(rx_file);
            % plot(real(raw_data)) check the bits toggled at ADC

        
        %% Reshape array into matrix of pulses
            pulse_matrix = reshape(raw_data,[length(raw_data)/number_pulses,number_pulses]); %reshape array to individual pulses
            adc_max = max(real(raw_data))
            adc_rms = rms(real(raw_data))
            dbfs = 10*log10(adc_rms*sqrt(2)/1)
            dbm = 30+ 20*log10(adc_rms/(sqrt(50)))

            clear raw_data
            figure
            t = linspace(0,0.5e-3,length(real(pulse_matrix(:,1))));
            plot(t,real(pulse_matrix(:,1)));
                ylabel('ADC value (V)')
                xlabel('Time (s)');      
                xlim([-inf +inf])
           pwr = abs(fftshift(fft(pulse_matrix(:,200),length(pulse_matrix(:,8)))));
%            pwr = pwelch(pulse_matrix(:,1),length(pulse_matrix(:,1)),0,'centered','power');
           f = linspace(-15,15,length(pulse_matrix(:,8)));
           figure
           plot(f,20*log10(pwr./max(pwr(:))))
           grid on; grid minor;
           xlabel('Frequency (MHz)');
           ylabel('Relative power (dB)');

                %                 title("Receive Channel Time Series - Max val: " + adc_max);


        %% Deramp Signal
            deramped_signal = zeros(size(pulse_matrix,1),size(pulse_matrix,2));
            tic
            deramped_signal = dechirp(pulse_matrix,refsig);
            clear pulse_matrix
            deramp_data = toc

        %% Decimate Signal 
        %%determine decimation factor number_pulses
        %decimation factor must be by an integer thus the maximum range will likley
        %be non integer due to rounding
            if_max = (max_range * 2 * slope)/C; % if_max is the maximum frequency were intrested in
            decimation_factor = Fs/if_max; % dec_fac is the non integer ideal decimation factor
            decimation_factor_rounded = floor(size(deramped_signal,1)/decimation_factor); %dec_fac_rnd is the closest integer 
            decimation_factor_actual = ceil(size(deramped_signal,1)/decimation_factor_rounded); %determines actual optmal decimation factor
            if_freq_actual = Fs/decimation_factor_actual; %determines subsequent maximum if frequency 
            actual_max_range = beat2range(if_freq_actual,slope); %translates if frequency to range
       % allocate memory 
            decimated_signal_size = size(decimate(deramped_signal(:,1),decimation_factor_actual),1);
            decimated_signal = zeros(decimated_signal_size,size(deramped_signal,2)); %initiate array

       % low-pass filter then decimate (requires alot of computation)
       if lp_filter == true
            tic
            for i=1:number_pulses
                 decimated_signal(:,i) = decimate(deramped_signal(:,i),decimation_factor_actual);
            end
            clear Data_Deramped
            Filter_and_Decimation_time = toc
        else 
            tic 
            for i=1:number_pulses
                 decimated_signal(:,i) = downsample(deramped_signal(:,i),decimation_factor_actual);
            end
            clear Data_Deramped
            Decimation_time = toc
        end
        
        
        return
    


end 




























if capture_duration > 62
    %If capture duration is over 30s, load capture in indiviual 30s sections to
    %process and concatanate at the end.
    window_size = 5; %window size in s
    number_windows = capture_duration/window_size;

        for window=0:number_windows

            %% load rx file
                rx_file= file_location + ".sc16q11";
                raw_data = load_sc16q11_window(rx_file,window,Fs,window_size);
                % plot(real(raw_data)) check the bits toggled at ADC

            %% Reshape array into matrix of pulses
                pulse_matrix = reshape(raw_data,[length(raw_data)/(number_pulses/number_windows),number_pulses/number_windows]); %reshape array to individual pulses
                clear raw_data
        %       figure 
        %       spectrogram(pulse_matrix(:,1000),128,100,100,Fs,'centered','yaxis')

            %% Deramp Signal
                deramped_signal = zeros(size(pulse_matrix,1),size(pulse_matrix,2));
                tic
                deramped_signal = dechirp(pulse_matrix,refsig);
                clear pulse_matrix
                deramp_data = toc

            %% Decimate Signal 

            if window == 0 
                %%determine decimation factor 
                %decimation factor must be by an integer thus the maximum range will likley
                %be non integer due to rounding
                    if_max = (max_range * 2 * slope)/C; % if_max is the maximum frequency were intrested in
                    decimation_factor = Fs/if_max; % dec_fac is the non integer ideal decimation factor
                    decimation_factor_rounded = round(size(deramped_signal,1)/decimation_factor); %dec_fac_rnd is the closest integer 
                    decimation_factor_actual = round(size(deramped_signal,1)/decimation_factor_rounded); %determines actual optmal decimation factor
                    if_freq_actual = Fs/decimation_factor_actual; %determines subsequent maximum if frequency 
                    actual_max_range = beat2range(if_freq_actual,Slope); %translates if frequency to range
            end
            %decimate signal
                decimated_signal = zeros(decimation_factor_rounded + 1,size(deramped_signal,2)); %initiate array
                tic
                for i=1:Num_pulses/number_windows
                decimated_signal(:,i) = decimate(deramped_signal(:,i),decimation_factor_actual);
                end
                clear Data_Deramped
                Decimation_time = toc

            %% FFT Signal
                tic
                processed_data = fft(decimated_signal,size(decimated_signal,1));
                clear Dec_Deramped 
            % if first winodow create empty matrix to store variables 
            window_width = size(processed_data,1);
            if window == 0
                Final_Data = zeros(number_windows*window_width,size(processed_data,2));
                Final_Data(1:window_width,:) = processed_data;
            else 
                Final_Data(window*window_width:(window+1)*window_width,:) = processed_data;   
            end
        end
end
    
end

