function [det_result,pd] = computePD(detections,ground_truth_range,ground_truth_sog,range_error_lim,sog_error_lim)
%COMPUTEPD compute the probability of detection using the ground truth data
%to compare to detections 
% detections: cell array where each cell contains the detections for each
% CPI. Column 1 is range, column 2 is Doppler velocity in m/s
number_cpi = size(detections,2);
det_result= zeros(number_cpi,1);
% loop through CPIs
for i=1:number_cpi
    target_range = ground_truth_range(i); % Evaluate ground truth range for CPI
    target_sog = ground_truth_sog(i);        % Evaluate ground truth SoG for CPI
    % loop through detections
    for j = 1:size(detections{i},2)
        % range error between ground truth and detection
        r_error = detections{i}(1,j) - target_range;
        % if range error is less than the range bin size
        if abs(r_error) < range_error_lim
            % is doppler less than the velocity bin size 
            d_error = detections{i}(2,j) - target_sog;
            if abs(d_error) < sog_error_lim
                det_result(i) = 1;
            end
        end
    end        
end
% Calculate the Pd
pd = (sum(det_result)/number_cpi)*100;
end