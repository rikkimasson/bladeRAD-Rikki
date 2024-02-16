function [signal, noise, SNR, slices_range_doppler, index_signal_max, noise_abs, slices_range_doppler_abs] = slice_SNR(range_doppler_slice, roi_indices_rows, roi_indices_cols)

%   This function takes a range doppler slice, a cell array containing
%   complex numbers, uses the ABS function in MATLAB to calculate the complex
%   numbers into absolute values, and then converts the absolute values into
%   signal power. Once this is known, for each slice, the Signal is identified
%   as the Peak power, and the Noise is calculated as the mean power (the average
%   of all the numbers in the set).  The SNR is then calculated from these.
%   All 3 are output.
%   A data masking has been hard-coded into this function. The first 5
%   columns of the data are set to zero in an attempt to mask the false
%   signal peaks observed at close range. 5 columns were chosen because
%   each column equates to a range bin, active.range_axis(5) = 19.9862 metres.
%   After viewing the data 20 metres seemed sensible.

%Create region of interest matrix
roi_data = range_doppler_slice(roi_indices_rows, roi_indices_cols);
complex_2_abs_roi = abs(roi_data);
complex_2_abs = abs(range_doppler_slice);
abs_2_power = 20*log10(complex_2_abs); % [dB]
abs_2_power_masked = abs_2_power;
abs_2_power_masked(:, 1:5) = 0; % Mask the first 5 columns = zero
[val, idx] = max(abs_2_power_masked,[],"all");
% [val, idx] = max(abs_2_power,[],"all");
signal = val; % [dB]
index_signal_max = idx;
noise = 20*log10(mean(complex_2_abs_roi(:))); % [dB]
noise_abs = mean(complex_2_abs_roi(:)); %  [linear]
%noise = 20*log10(mean(complex_2_abs(:)));
SNR = signal - noise; % [dB]
slices_range_doppler = abs_2_power; % [dB]
slices_range_doppler_abs = complex_2_abs; % [linear]
end

% function [signal, noise, SNR, slices_range_doppler, index_signal_max] = slice_SNR(range_doppler_slice)
% 
% %   This function takes a range doppler slice, a cell array containing
% %   complex numbers, uses the ABS function in MATLAB to calculate the complex
% %   numbers into absolute values, and then converts the absolute values into
% %   signal power. Once this is known, for each slice, the Signal is identified
% %   as the Peak power, and the Noise is calculated as the mean power (the average
% %   of all the numbers in the set).  The SNR is then calculated from these.
% %   All 3 are output.
% 
% complex_2_abs = abs(range_doppler_slice);
% abs_2_power = 20*log10(complex_2_abs);
% %signal = max(abs_2_power(:));
% [val, idx] = max(abs_2_power(:));
% signal = val;
% index_signal_max = idx;
% noise = 20*log10(mean(complex_2_abs(:)));
% SNR = signal - noise;
% slices_range_doppler = abs_2_power;
% end