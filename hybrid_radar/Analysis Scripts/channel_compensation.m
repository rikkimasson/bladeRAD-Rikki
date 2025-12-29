function [XF_comp,Comps] = channel_compensation(XF,CS, carrier_locations, numb_carriers)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    XF_comp=XF;
    Comps=complex(zeros(1,length(XF_comp)));
    for i=1:1:numb_carriers
        % if i==3408
        %     jsfdhf=3;
        % end
        [~,II]=min(abs(i-carrier_locations));
        if II==351
          II=352;
        end
        XF_comp(i)=XF(i)./CS(II);
        Comps(i)=CS(II);

    end
end