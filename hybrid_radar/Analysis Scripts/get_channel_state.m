function [CS] = get_channel_state(X_pilots_int,prbs_seq,carrier_locations)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
pilot_locations=(4/3)*2*(0.5-prbs_seq)+1j*0;

CS=X_pilots_int./pilot_locations(carrier_locations)';
end