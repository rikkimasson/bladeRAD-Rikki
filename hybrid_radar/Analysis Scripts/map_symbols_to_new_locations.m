function [perfect_symbols] = map_symbols_to_new_locations(i_symbol,q_symbol,numb_carriers,prbs_seq,sypo,carrier_locations,tps_carriers)
%UNTITLED10 Summary of this function goes here
perfect_symbols=complex(zeros(1,numb_carriers));
for i=1:1:numb_carriers
    if any( carrier_locations== i) 
        % continue 
        perfect_symbols(i)=(4/3)*2*(0.5-prbs_seq(i))+1j*0;
    elseif any(tps_carriers==i)
        if sypo(i_symbol(i))<0
        perfect_symbols(i)=-1;
        else 
            perfect_symbols(i)=1;
        end
    else
        perfect_symbols(i)=sypo(i_symbol(i))+1j*sypo(q_symbol(i));
    end
end
end