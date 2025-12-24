function [i_symbol,q_symbol] = get_actual_symbols(XF,num_carriers,sypo,carrier_locations)

i_symbol=zeros(1,num_carriers);
q_symbol=zeros(1,num_carriers);

for i=1:1:num_carriers
    if any( carrier_locations== i)
        continue 
    else
        [~,i_symbol(i)]= min(abs(real(XF(i))-sypo));
        [~,q_symbol(i)]= min(abs(imag(XF(i))-sypo));
    end
end
end