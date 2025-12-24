function [output_signal] = produce_ideal_ofdm_symbol(input_signal,S,D,numb_carriers,permanant_carriers)
    % takes in reference (or surveillance signal and demodulates it and
    % recreates the signal based on the actual bits, produces a cleaned
    % signal which can be used as reference
    output_signal=complex(zeros(1,length(input_signal)));
    offset_order=[10,13,16,19];
    offset_spacing=12;
    offset_probability=0;
    guess_offset=1;

    prbs_seq = generate_prbs(numb_carriers);

    scale_factor=6.3636;
    sypo=[-7,-5,-3,-1,1,3,5,7]/scale_factor;

    symbol_starts=find_symbol_starts(input_signal,S,D);

    for i=1:1:length(symbol_starts)
        msignal=input_signal(symbol_starts(i)+D+1:symbol_starts(i)+D+S+1);
        XX=linspace(-1-5+7,6819-5+4,length(msignal));
        XV=linspace(1,numb_carriers,numb_carriers);
        new_signal=interp1(XX,msignal,XV,"spline");


        XF_int=fft(new_signal,numb_carriers);

        offset=get_carrier_offset(XF_int(1:60),guess_offset,offset_probability,offset_order,offset_spacing);
        
        guess_offset=update_offset(offset);

        carrier_locations=get_all_carrier_locations(offset,permanant_carriers,numb_carriers);
        
        X_pilots_int=XF_int(carrier_locations);
        [CS] = get_channel_state(X_pilots_int,prbs_seq,carrier_locations);
        
        [XF_comp] = channel_compensation(XF_int,CS, carrier_locations, numb_carriers);
            
        [i_symbol,q_symbol] = get_actual_symbols(XF_comp,numb_carriers,sypo,carrier_locations);
        
        [XF_ideal] = map_symbols_to_new_locations(i_symbol,q_symbol,numb_carriers,prbs_seq,sypo,carrier_locations);
        figure
        hold on
        scatter(real(XF_comp(1:1500)),imag(XF_comp(1:1500)))
        
        %should compensate for the frequency offset which I can get from
        %the slope of the plot of the carrier symbols
        [output_IQ] = produce_IQ_data(XF_ideal,numb_carriers,length(msignal));
        
        % check that this is in fact aligned correctly and does not need a
        % +1 or -1 anywhere 
        output_signal(symbol_starts(i)+D+1:symbol_starts(i)+D+S+1)=output_IQ;
        output_signal(symbol_starts(i):symbol_starts(i)+D+1)=output_IQ(S-D:end);
    end




end