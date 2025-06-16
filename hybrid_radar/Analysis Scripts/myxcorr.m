function [out] = myxcorr(input1,input2)

l1=length(input1);
l2=length(input2);
l3=l1+l2;
out=(zeros(1,l1+l2));

extend1=complex(zeros(1,l3));
extend2=complex(zeros(1,l3));
extend1(1:l1)=input1;
extend2(l3-l2+1:l3)=input2;

for i=1:1:l3
    out(i)=abs(sum(circshift(extend1,i).*conj(extend2)));
    % if i==22987
        % figure
        % plot(abs(circshift(extend1,i).*conj(extend2)))
    % end
end

figure,plot(out)

figure,
hold on
plot(real(out))
plot(imag(out))

end