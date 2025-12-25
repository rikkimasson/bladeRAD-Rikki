function [symbols_starts,symbol_phases] = find_symbol_starts(input,S,D,dt)
% looks for start points of the ofdm symbols
P1=(zeros(1,length(input)));
    for i=1:1:length(input)-S-D-1
        P1(i)=(sum(conj(input(i:i+D)).*input(i+S:i+S+D)));
    end
    
    % figure
    % plot(abs(P1))
    % 
    % figure
    % plot(angle(P1))
    P1_angle=angle(P1);
    P1=abs(P1);

    % t=linspace(0,dt*length(P1),length(P1));
    % corr=exp(1j*0.150455/S/dt*t)
    % 
    % temp=input.*corr;
    % 
    % for i=1:1:length(input)-S-D-1
    %     P1_temp(i)=(sum(conj(temp(i:i+D)).*temp(i+S:i+S+D)));
    % end

    % figure
    % plot(angle(P1_temp))
    

    % 
    symbols_starts=[];
    % symbol_phases=[]
    [m,I_new]=max(P1(1:2*S));
    symbols_starts=[I_new];
    while 1
        if I_new+D+S+200> length(P1)
            break
        end
        I_prev=I_new;
        [m_new, I_new]=max(P1(I_prev+D+S-1000:I_prev+D+S+1000));
        I_new=I_new+I_prev+D+S-1001;
        if abs(m_new-m)/m<0.5 
            symbols_starts=[symbols_starts,I_new];
        else 
            fprintf('something has gone wrong with tracking begining of frame')
        end

    end
    
    symbol_phases=zeros(1,length(symbols_starts));
    for i=1:1:length(symbols_starts)
        symbol_phases(i)=P1_angle(symbols_starts(i));
    end

    % state=0;
    % mthreshold=0.75*max((P1));
    % down_threshold=0.2*max((P1));
    % symbols_starts=[];
    % for i=1:1:length(P1)
    %     if (P1(i))>mthreshold
    %         if state==0
    %             current_peak=(P1(i));
    %             current_index=i;
    % 
    %         elseif abs(P1(i))>current_peak
    %             current_peak=(P1(i));
    %             current_index=i;
    %         end
    %         state=1;
    % 
    %     else
    %         if (P1(i))<down_threshold && state==1
    %             symbols_starts=[symbols_starts,current_index];
    %              state=0;
    %             current_peak=0;
    %         end
    %         % if state==1
    %         %     symbols_starts=[symbols_starts,current_index];
    %         % end
    %         % state=0;
    %         % current_peak=0;
    %     end
    % end




end