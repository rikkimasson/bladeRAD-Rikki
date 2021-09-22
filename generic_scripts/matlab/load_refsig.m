function refsig = load_refsig(BW,Fc,pulse_dur)
%LOAD_REFSIG load reference signal to match filter or deramp signal
pulse_dur = pulse_dur*1e3;
Fc = Fc*1e-9
file_name = string(pulse_dur) + 'ms_' + string(Fc) + 'GHz_' + BW + 'MHz.sc16q11';
refsig = load_sc16q11(file_name);

end
