function [x3] = tone_generator(fDS,Pulse_dur,Fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Samp_dur = 1/Fs;

t = 0:Samp_dur:Pulse_dur; %time vector of chirp
x3 = exp(1i*2*pi*(fDS)*t(1:length(t))); %up chirp from 0 to BW/2



end

