function [x3] = saw_LFM_chirp(Bw,Pulse_dur,Fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Samp_dur = 1/Fs;

t = 0:Samp_dur:Pulse_dur; %time vector of chirp
NUM_SAMPLES = length(t)-1;
x1 = exp(1i*pi*(Bw/Pulse_dur)*t(1:length(t)/2).^2); %up chirp from 0 to BW/2
x2 =fliplr(x1); %up chirp from -BW/2 to 0 
x3 = 0.9 * [x2 x1]; %concatenate two two vectors



end

