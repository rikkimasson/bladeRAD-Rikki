function detection = CFAR(matched_data,n_pulses,n_rbins,pulse_no)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%CFAR Design parametes
        threshold_offset = 3;
        nguard_cells = 10;  %number of guard cells at centre of window
        nref_cells = 20;    %number of refence cells at edge of window
        window_l = 2 *(nref_cells + nguard_cells ) + 1; %total window length

%create vector to window the data with 
        window = zeros(1,window_l);
        window(1,1:nref_cells) = 1;
        window(1,2 + nref_cells + nguard_cells:end) = 1;

%% Create Cell Averages 

threshold_matrix = zeros(size(matched_data));       
        
for i = 1:n_pulses  % stride pulses
    
    dwell = abs(matched_data(:,i));
       
    for j = 1:n_rbins % stride range bins
                
        if j < nref_cells + nguard_cells % at start of dwell 
            
            overlap = nguard_cells+nref_cells - j; 
            crop_window_l = window_l - overlap;
%           j = j
            threshold_matrix(j,i) =  dot( dwell(1+j:1+j+crop_window_l) , window(overlap:end));
                if overlap+1 > nref_cells
                    threshold_matrix(j,i) = threshold_matrix(j,i) / nref_cells;
                else
                    threshold_matrix(j,i) = threshold_matrix(j,i) / (nref_cells - overlap);
                end
            
        elseif j > n_rbins - (nref_cells + nguard_cells+1)  % at end of dwell
            
%             j = j 
            overlap = nguard_cells+nref_cells - (n_rbins - j); 
            crop_window_l = window_l - overlap;
            dwell_p = dwell(j-nref_cells-nguard_cells:j+nref_cells + nguard_cells-overlap);
            dwell_pl = length(dwell_p);
            threshold_matrix(j,i) =  dot( dwell(j-nref_cells-nguard_cells:j+nref_cells + nguard_cells-overlap) , window(1:window_l-overlap));
                if overlap+1 > nref_cells
                    threshold_matrix(j,i) = threshold_matrix(j,i) / nref_cells;
                else
                    threshold_matrix(j,i) = threshold_matrix(j,i) / (nref_cells - overlap);
                end
            
        else % all full window lengths
             
             dwell_p = dwell(1+j-nref_cells-nguard_cells:1+j+nref_cells+nguard_cells);
             dwell_pl = length(dwell_p);
%              j = j
             threshold_matrix(j,i) =  dot(dwell((1+j-(nref_cells + nguard_cells)):1+j+(nref_cells + nguard_cells)),window); % multiply and sum window by dwell period under evaluation
             threshold_matrix(j,i) = threshold_matrix(j,i)  / (2 * nref_cells); %average sum by the number of refernce cells 
             
        end
    end   
end    



%% plot indvidual pulses
x = transpose(linspace(1,n_rbins,n_rbins));
plot(x,20*log10(abs(matched_data(:,pulse_no))),x,20*log10(threshold_matrix(:,pulse_no)));

%% Target decsion 
threshold_matrix = threshold_matrix * threshold_offset; % Increase threshold by integer
decsion_matrix = zeros(size(matched_data));       
x = transpose(linspace(1,n_rbins,n_rbins));
plot(x,20*log10(abs(matched_data(:,pulse_no))),x,20*log10(threshold_matrix(:,pulse_no)));

for i = 1:n_pulses
    for j = 1:n_rbins
        if threshold_matrix(j,i) < matched_data(j,i)
            decsion_matrix(j,i) = 1*1000;
        else 
            decsion_matrix(j,i) = 0;
        end
    end
end
detection = decsion_matrix;
end

            