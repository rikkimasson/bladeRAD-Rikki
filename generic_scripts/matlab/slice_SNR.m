function [signal, noise, SNR, slice_range_doppler, index_signal_max] = slice_SNR(range_doppler_slice)

%   This function takes a range doppler slice, a cell array containing
%   complex numbers, uses the ABS function in MATLAB to calculate the complex
%   numbers into absolute values, and then converts the absolute values into
%   signal power. Once this is known, for each slice, the Signal is identified
%   as the Peak power, and the Noise is calculated as the mean power (the average
%   of all the numbers in the set).  The SNR is then calculated from these.
%   All 3 are output.

complex_2_abs = abs(range_doppler_slice);
abs_2_power = 20*log10(complex_2_abs);
%signal = max(abs_2_power(:));
[val, idx] = max(abs_2_power(:));
signal = val;
index_signal_max = idx;
noise = 20*log10(mean(complex_2_abs(:)));
SNR = signal - noise;
slice_range_doppler = abs_2_power;
end