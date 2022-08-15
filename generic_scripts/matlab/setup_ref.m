function [exit_code] = setup_ref(SDR_No,clock_ref)
%create_shell_command - Used to create commands to run shell scripts for
%each SDR.

sdr_serial = select_SDR(SDR_No); % find SDR serial number

% run async in terminal window or run in matlab
    command = "~/repos/bladeRAD/generic_scripts/initialise_clk_dist.sh " + sdr_serial + " " + clock_ref; 
    exit_code = system(command);
