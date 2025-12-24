function [all_carriers] = get_all_carrier_locations(offset, get_all_carrier_locations,numb_carriers)
%UNTITLED6 Summary of this function goes here
    if offset==1
        start=10;
    elseif offset==2
        start=13;
    elseif offset==3
        start=16;
    elseif offset ==4
        start=19;
    else
        start=10;
        printf('this should not happen')
    end
    temp=start:12:numb_carriers;

    all_carriers=unique([temp,get_all_carrier_locations]);

end