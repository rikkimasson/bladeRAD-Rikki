function refsig = load_refsig(Bw,Fs,Fc,pulse_duration)
%LOAD_REFSIG load reference signal to match filter or deramp signal
cli_pulse_duration = pulse_duration*1e3;
Fc = Fc*1e-9;
file_name = string(cli_pulse_duration) + 'ms_' + string(Fc) + 'GHz_' + Bw + 'MHz.sc16q11';
try
    refsig = load_sc16q11(file_name);
catch
    'No pre-recorded refsig'
    Bw = Bw *1e6;
    refsig = transpose(saw_LFM_chirp(Bw,pulse_duration,Fs));
end
