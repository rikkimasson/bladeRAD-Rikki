function [trig_flag,command] = create_shell_command(popout, test_id, cap_samps, pulses, delay, tx_gain, rx1_gain, rx2_gain, center_freq, bw, SDR_No, trigger,clock,tx_rx)
%create_shell_command - Used to create commands to run shell scripts for
%each SDR.

if trigger == "master"
    trig_flag = true;
else 
    trig_flag = false;
end

sdr_serial = select_SDR(SDR_No); % find SDR serial number

% run async in terminal window or run in matlab
if popout == true
    command = "gnome-terminal -- ~/repos/bladeRAD/generic_scripts/bladeRAD_shell_script.sh " + test_id + " " + cap_samps + " " + pulses + " " + delay + " " + tx_gain + " " + rx1_gain + " " + rx2_gain + " " + center_freq + " " + bw + " " + sdr_serial + " " + trigger + " " + clock + " " + tx_rx; 
else
    command = "~/repos/bladeRAD/generic_scripts/bladeRAD_shell_script.sh " + test_id + " " + cap_samps + " " + pulses + " " + delay + " " + tx_gain + " " + rx1_gain + " " + rx2_gain + " " + center_freq + " " + bw + " " + sdr_serial + " " + trigger + " " + clock + " " + tx_rx; 

end

end
