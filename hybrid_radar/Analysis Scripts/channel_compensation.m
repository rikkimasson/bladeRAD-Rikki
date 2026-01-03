function [XF_comp,Comps] = channel_compensation(XF,CS, carrier_locations, numb_carriers)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
figure
plot(abs(CS))

figure
plot(angle(CS))

temp=abs(CS);
dont_use_carrier=[];

% temp2=zeros(1,length(CS));

for i=2:1:length(CS)-1
    fl_avg=0.5*(temp(i-1)+temp(i+1));
    % temp2(i)=abs(temp(i)-fl_avg)/fl_avg;
    if abs(temp(i)-fl_avg)/fl_avg>0.1
        % abs(temp(i)-fl_avg)/fl_avg
        dont_use_carrier=[dont_use_carrier,i];
    end
end
% 
% figure
% plot(temp2)
good_carriers=carrier_locations;
CS_good=CS;
for i=1:1:length(dont_use_carrier)
    good_carriers(dont_use_carrier(i))=[];
    CS_good(dont_use_carrier(i))=[];
end

    XF_comp=XF;
    Comps=complex(zeros(1,length(XF_comp)));
    for i=1:1:numb_carriers
        if any( carrier_locations== i)
            [~,II]=min(abs(i-carrier_locations));
            XF_comp(i)=XF(i)./CS(II);
            Comps(i)=CS(II);

        else
            % [~,II_l]=max((i-good_carriers));
            % [~,II_h]=max((good_carriers-i));
             [lower, upper] = findClosestNeighbors(good_carriers, i);
            CS_interp=interp1([good_carriers(lower),good_carriers(upper)],[CS_good(lower),CS_good(upper)],i);
            XF_comp(i)=XF(i)./CS_interp;
            Comps(i)=CS_interp;

        end
        % if i==3408
        %     jsfdhf=3;
        % end
        % [~,II]=min(abs(i-carrier_locations));
        % if II==351
        %   II=352;
        % end
        % XF_comp(i)=XF(i)./CS(II);
        % Comps(i)=CS(II);

    end
end