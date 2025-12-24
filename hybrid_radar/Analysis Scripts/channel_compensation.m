function [XF_comp] = channel_compensation(XF,CS, carrier_locations, numb_carriers)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    XF_comp=XF;
    for i=1:1:numb_carriers
        [~,II]=min(abs(i-carrier_locations));
        XF_comp(i)=XF(i)./CS(II);

    end
end