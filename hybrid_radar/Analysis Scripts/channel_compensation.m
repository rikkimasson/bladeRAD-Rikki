function [XF_comp] = channel_compensation(XF,CS, carrier_locations, numb_carriers)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    for i=1:1:6817
        [M,II]=min(abs(i-carrier_locations));
        XF_comp(i)=XF_int(i)./CS(II);

    end
end