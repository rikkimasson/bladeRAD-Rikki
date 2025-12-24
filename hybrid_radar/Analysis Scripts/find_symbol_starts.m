function [symbols_starts] = find_symbol_starts(input,S,D)
% looks for start points of the ofdm symbols
P1=(zeros(1,length(input)));
    for i=1:1:length(input)-S-D-1
        P1(i)=abs(sum(conj(input(i:i+D)).*input(i+S:i+S+D)));
    end
    

    % 
    symbols_starts=[];
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