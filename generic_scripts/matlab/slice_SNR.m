function [signal, noise, SNR] = slice_SNR(range_doppler_slice)

%   This function takes a range doppler slice, a cell array contains
%   complex numbers, uses the ABS function in MATLAB to calculate the complex
%   numbers into absolute values, and then converts the absoute values into
%   signal power. Once this is known, for each slice, the Signal is identified
%   as the Peak power, and the noise is calculated as the mean power (the average
%   of all the numbers in the set).  The SNR is then calculated from these.
%   All 3 are output.

complex_2_abs = abs(range_doppler_slice);
abs_2_power = 20*log10(complex_2_abs);
signal = max(abs_2_power(:))
noise = 20*log10(mean(complex_2_abs(:)))
SNR = signal - noise

end