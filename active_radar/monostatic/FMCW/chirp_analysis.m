          
        [spect,f] = spectrogram(refsig,127,100,127,active.Fs,'centered','yaxis');
        time_axis = linspace(0,active.pulse_duration,size(refsig,2));
        spect= 10*log10(abs(spect./max(spect(:))));
        figure
%         subplot(1,2,2)
        imagesc(time_axis,f*10^-6,(spect),[-50 0]);
        ax = gca;
        ax.YDir = 'normal'
        ylim([-15 15])
        xlabel('Time (\mus)')
        ylabel('Frequency (MHz)')  
        % ylabel('Radial Velocity (mph)')   
        ylabel('Frequency (MHz)')  
        fig_title = "Monostatic Spectrogram - Test " + Experiment_ID;
        c = colorbar
        c.Label.String='Normalised Power (dB)';
        
        
        figure
%         subplot(1,2,1)
        plot((time_axis),real(refsig));
        grid on 
        grid minor
        xlabel('Time (\mus)')
        ylabel('Amplitude')  
%         xlim([0 200])

        freq = abs(fftshift(fft(refsig,size(refsig,2)*2)));
        freq = 10*log10(abs(freq./max(freq(:))));
        f = linspace(-15,15,size(freq,2));
        subplot(1,2,2)
        plot(f,freq);
%         xlim([-16 16])
%         ylim([-20 5])
        xlabel('Frequency (MHz)')
        ylabel('Normalised Power (dB)')  
        grid on
        grid minor