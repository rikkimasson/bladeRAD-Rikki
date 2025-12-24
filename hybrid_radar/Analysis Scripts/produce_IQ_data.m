function [res2_signal_per] = produce_IQ_data(perfect_symbols,numb_carriers,sample_length)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
temp4_per=ifftshift(perfect_symbols);
temp5_per=(sample_length/numb_carriers)*[temp4_per(1:3409),zeros(1,sample_length-numb_carriers),temp4_per(3410:end)];
res2_signal_per=ifft(temp5_per);
end