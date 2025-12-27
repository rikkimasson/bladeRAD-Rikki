function [output_offset] = get_carrier_offset(XF_int,guess_offset, offset_probability,offset_order,offset_spacing)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    % [B,I]=sort(abs(XF_int));

    XF=abs(XF_int);
    temp=zeros(1,4);
    temp(1)=XF(10)+XF(offset_order(1)+offset_spacing)+XF(offset_order(1)+2*offset_spacing)+XF(offset_order(1)+3*offset_spacing);
    temp(2)=XF(13)+XF(offset_order(2)+offset_spacing)+XF(offset_order(2)+2*offset_spacing)+XF(offset_order(2)+3*offset_spacing);
    temp(3)=XF(16)+XF(offset_order(3)+offset_spacing)+XF(offset_order(3)+2*offset_spacing)+XF(offset_order(3)+3*offset_spacing);
    temp(4)=XF(19)+XF(offset_order(4)+offset_spacing)+XF(offset_order(4)+2*offset_spacing)+XF(offset_order(4)+3*offset_spacing);

    [B,I]=max(temp);

    temp2=temp;
    temp2(I)=[];

    second=max(temp2);

    myguess=(abs(B-second)/second)*100;

    if (guess_offset==I)
        output_offset=I;
    elseif (offset_probability>myguess)
        output_offset=guess_offset;
        fprintf("metrics disagree on what is carrier offset\n")
    else
        output_offset=I;
        fprintf("metrics disagree on what is carrier offset\n")
    end

end