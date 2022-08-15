#!/bin/sh

sdr_serial=${1}
clock=${2}

IFS=$''


#create clock ref message
if [ "$clock" = 1 ]; then clock_ref='set clock_sel external';
elif [ "$clock" = 2 ]; then clock_ref='set clock_ref enable';
else clock_ref=''; fi




bladeRF-cli -d "*:serial=$sdr_serial" -e ' 	
    			'$clock_ref';
    			set clock_out enable;
			print'

