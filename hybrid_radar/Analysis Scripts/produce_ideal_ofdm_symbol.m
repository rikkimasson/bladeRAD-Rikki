function [output_signal,total_offsets] = produce_ideal_ofdm_symbol(input_signal,S,D,numb_carriers,permanant_carriers,tps_carriers,dt)
    % takes in reference (or surveillance signal and demodulates it and
    % recreates the signal based on the actual bits, produces a cleaned
    % signal which can be used as reference
    output_signal=complex(zeros(1,length(input_signal)));
    offset_order=[10,13,16,19];
    offset_spacing=12;
    offset_probability=0;
    offset=4;

    prbs_seq = generate_prbs(numb_carriers);

    scale_factor=6.3636;
    sypo=[-7,-5,-3,-1,1,3,5,7]/scale_factor;

    [symbol_starts,symbol_phases]=find_symbol_starts(input_signal,S,D,dt);
    save("symbol_starts_data_1","symbol_starts","symbol_phases")

    % load("symbol_starts_data")
    %


    % temp_diff=symbol_starts(2:length(symbol_starts))-symbol_starts(1:length(symbol_starts)-1);    
    % figure
    % plot(temp_diff)
    %looks like the packet is streaming continuously so can update offset
    %automatically


    % figure
    % hold on
    total_offsets=zeros(1,length(symbol_starts));

    for i=1:1:length(symbol_starts)
        msignal=input_signal(symbol_starts(i)+D+1:symbol_starts(i)+D+S+1);
        check_me=msignal;
        t=linspace(0,dt*length(msignal),length(msignal));
        corr=exp(-1j*symbol_phases(i)/S/dt*t);
        msignal=msignal.*corr;
        

       
        XX=linspace(-1-5+7,6819-5+4,length(msignal));
        XV=linspace(1,numb_carriers,numb_carriers);
        new_signal=interp1(XX,msignal,XV,"spline");


        XF_int=fft(new_signal,numb_carriers);

        XF_int=fftshift(XF_int);

        % figure 
        % hold on
        % scatter(real(XF_int(1:200)),imag(XF_int(1:200)))
        guess_offset=offset;
        guess_offset=update_offset(guess_offset);
        if i>1
             if abs((abs(symbol_starts(i)-symbol_starts(i-1))/13440)-1)<0.2
                offset_probability=100;
             else
                 offset_probability=0;
                 fprintf("this should not really happen\n");
             end
        end

        offset=get_carrier_offset(XF_int(1:60),guess_offset,offset_probability,offset_order,offset_spacing);
        
        total_offsets(i)=offset;

        carrier_locations=get_all_carrier_locations(offset,permanant_carriers,numb_carriers);
        
        X_pilots_int=XF_int(carrier_locations);

        % figure
        % hold on
        % plot(unwrap(angle(X_pilots_int)))
        % 
        % figure
        % plot((angle(X_pilots_int)))
        
        % seems like I should probably use a linear interpolation of the
        % pilot symbols because otherwise will get the intersymbol
        % interference material problems
        [CS] = get_channel_state(X_pilots_int,prbs_seq,carrier_locations);
        
        [XF_comp,Comps] = channel_compensation(XF_int,CS, carrier_locations, numb_carriers);
            
        [i_symbol,q_symbol] = get_actual_symbols(XF_comp,numb_carriers,sypo,carrier_locations,tps_carriers);
        
        [XF_ideal] = map_symbols_to_new_locations(i_symbol,q_symbol,numb_carriers,prbs_seq,sypo,carrier_locations,tps_carriers);
        % figure
        % hold on
        % scatter(real(XF_int),imag(XF_int))
        % % scatter(real(XF_int(indexess)),imag(XF_int(indexess)),'filled')
        
        % figure,
        % hold on
        % scatter(real(XF_ideal),imag(XF_ideal))


        % figure,
        % hold on
        % scatter(real(XF_comp),imag(XF_comp))
        % scatter(real(XF_comp(3405:3410)),imag(XF_comp(3405:3410)),'filled')
        % scatter(real(XF_comp(3405)),imag(XF_comp(3405)),'filled')
        % scatter(real(XF_comp(3410)),imag(XF_comp(3410)),'filled')
        % scatter(real(XF_comp(tps_carriers)),imag(XF_comp(tps_carriers)),'filled')
        % scatter(real(XF_comp(3410:3413)),imag(XF_comp(3410:3413)),'filled')
        % scatter(real(XF_comp(indexess)),imag(XF_comp(indexess)),'filled')
        % x_check=[-1.21,-1.289,-.557,-0.2616];
        % y_check=[0.255,-0.788,-0.466,-1.223];
        % x_check=0.057;%-0.0638;%-0.3923;%-0.292;
        % y_check=-0.266;%-0.078;%-0.039;%-0.53;
        % epsilon=0.0112;
        % indx=(abs(real(XF_comp)-x_check(4))<epsilon & abs(imag(XF_comp)-y_check(4)));
        % find(indx)
        % sum(indx)

        % indexess=0.05>imag(XF_comp) & imag(XF_comp)>-0.075 & real(XF_comp)>-1.05 & real(XF_comp)<-0.85
        % sum(indexess)
        % indexess
        % figure,scatter(real(XF_comp(indexess)),imag(XF_comp(indexess)))

        % figure,plot(indexess)
        %should compensate for the frequency offset which I can get from
        %the slope of the plot of the carrier symbols
        % XF_ideal=XF_ideal.*Comps;
        % [output_IQ] = produce_IQ_data(XF_ideal.*Comps,numb_carriers,length(msignal),CS);
        % 
        XF_comp=XF_comp.*Comps;
        [output_IQ] = produce_IQ_data(XF_comp,numb_carriers,length(msignal),CS);

        corr=exp(1j*symbol_phases(i)/S/dt*t);
        output_IQ=output_IQ.*corr;

        % figure
        % hold on
        % plot(real(check_me)/max(abs(real(check_me))))
        % plot(real(output_IQ)/max(abs(real(output_IQ))))
        % 
        % figure
        % hold on
        % plot(unwrap(angle(check_me)))
        % plot(unwrap(angle(output_IQ)))
        % 

        % temp4_per=ifftshift(XF_int);
        % temp5_per=(length(msignal)/numb_carriers)*[temp4_per(1:3409),zeros(1,length(msignal)-numb_carriers),temp4_per(3410:end)];
        % XF_int_IQ=ifft(temp5_per);
        % 
        % XF_redoCS=XF_comp.*Comps;
        % 
        % temp4_per=ifftshift(XF_redoCS);
        % temp5_per=(length(msignal)/numb_carriers)*[temp4_per(1:3409),zeros(1,length(msignal)-numb_carriers),temp4_per(3410:end)];
        % XF_redoCS_IQ=ifft(temp5_per);
        % 
        % XF_ideal_CS=XF_ideal.*Comps;
        % 
        % temp4_per=ifftshift(XF_ideal_CS);
        % temp5_per=(length(msignal)/numb_carriers)*[temp4_per(1:3409),zeros(1,length(msignal)-numb_carriers),temp4_per(3410:end)];
        % XF_ideal_CS_IQ=ifft(temp5_per);


        % figure
        % hold on
        % plot(real(output_IQ)/max(abs(real(output_IQ))))
        % plot(real(msignal)/max(abs(real(msignal))))
        % figure
        % hold on
        % plot(real(check_me)/max(abs(real(check_me))))
        % plot(real(msignal)/max(abs(real(msignal))))
        % % plot(real(output_IQ)/max(abs(real(output_IQ))))
        % % plot(real(output_IQ_ideal)/max(abs(real(output_IQ_ideal))))
        % % plot(real(XF_int_IQ)/max(abs(real(XF_int_IQ))))
        % % plot(real(XF_redoCS_IQ)/max(abs(real(XF_redoCS_IQ))))
        % plot(real(XF_ideal_CS_IQ)/max(abs(real(XF_ideal_CS_IQ))))

        % plot(real(check_this)/max(abs(real(check_this))))
        % 
        % temp=xcorr(output_IQ,msignal);
        % 
        % figure,
        % plot(abs(temp))
       

        mmm=max(abs(check_me));
        m_small=max(abs(output_IQ));
        output_IQ=mmm/m_small*output_IQ;

        % figure
        % hold on
        % plot(real(check_me))
        % plot(real(output_IQ))
        
        output_signal(symbol_starts(i)+D+1:symbol_starts(i)+D+S+1)=output_IQ;
        output_signal(symbol_starts(i)-1:symbol_starts(i)+D)=output_IQ(S-D:end);

        % figure,
        % hold on
        % plot(real(output_signal(symbol_starts(i)-1:symbol_starts(i)+D+S+1)))
        % plot(real(input_signal(symbol_starts(i)-1:symbol_starts(i)+D+S+1)))

        % output_signal(symbol_starts(i)+D+1:symbol_starts(i)+D+S+1)=msignal;
        % output_signal(symbol_starts(i)-1:symbol_starts(i)+D)=msignal(S-D:end);
        % plot(real(output_signal))
    end




end