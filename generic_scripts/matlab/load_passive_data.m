function [ref_channel,sur_channel] = load_passive_data(file_location)
%PROCESS_PASSIVE_DATA Summary of this function goes here
%   Detailed explanation goes here


%% Load Signal    
    %load rx file
        rx_file = file_location + ".sc16q11";
        raw_data= load_sc16q11(rx_file); 
    
    % De-interleave data 
        ref_channel  = raw_data(1:2:end,:);
        sur_channel  = raw_data(2:2:end,:);
        clear raw_data

end

