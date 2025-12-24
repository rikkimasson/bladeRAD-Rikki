function [outputArg1,outputArg2] = produce_ideal_ofdm_symbol(input_signal,S,D,numb_carriers,permanant_carriers)
    % takes in reference (or surveillance signal and demodulates it and
    % recreates the signal based on the actual bits, produces a cleaned
    % signal which can be used as reference
    offset_order=[10,13,16,19];
    offset_spacing=12;
    offset_probability=0;
    guess_offset=1;

    prbs_seq = generate_prbs(numb_carriers);

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
        
        [XF_comp] = channel_compensation(XF_int,CS, carrier_locations, numb_carriers)


    end




end