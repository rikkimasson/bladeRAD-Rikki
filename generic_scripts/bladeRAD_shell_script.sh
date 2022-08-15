#!/bin/sh

test_id=$1
cap_samps=$2
pulses=$3
delay=$4
tx_gain=$5
rx1_gain=$6
rx2_gain=$7
center_freq=$8
bw=$9
sdr_serial=${10}
trigger=${11}
clock=${12}
tx_rx=${13}


IFS=$''

#Create trigger message
if [ "$trigger" = 'master' ]
then 
	triggerctrl="j51-1";
	fire="fire"
	chain=$tx_rx
else 	
	triggerctrl=""
	fire=""
	chain=""
fi	

export triggerctrl

#creatr clock ref message
if [ "$clock" = 1 ]; then clock_ref='set clock_sel external';
elif [ "$clock" = 2 ]; then clock_ref='set clock_ref enable';
else clock_ref=''; fi



# if Transmission ONLY is required run the following
if [ "$tx_rx" = 'tx' ]
then 
	bladeRF-cli -d "*:serial=$sdr_serial" -v verbose -e ' 
			
			set frequency tx '$center_freq'M;
			set samplerate tx '$bw'M;
			set bandwidth tx '$bw'M;
			set gain tx1 '$tx_gain';

			tx config file=/tmp/chirp.sc16q11 format=bin samples=16384 buffers=10 xfers=5 repeat='$pulses' timeout=60s delay='$delay' ;
	
    			'$clock_ref';
    			set clock_out enable;
			trigger j51-1 tx '$trigger';
			print;
            
            
			tx start;
			trigger '$triggerctrl' '$chain' '$fire';
			tx wait'
fi



# if Reception on Rx1 channel is required run the following
if [ "$tx_rx" = 'rx' ]
then 
	bladeRF-cli -d "*:serial=$sdr_serial" -v verbose -e '
			
			set frequency rx '$center_freq'M;
			set samplerate rx '$bw'M;
			set bandwidth rx '$bw'M;
			set agc off; 
			set gain rx1 '$rx1_gain' ; 

			
			rx config file=/tmp/active_'$test_id'.sc16q11 format=bin n='$cap_samps' samples=134144 buffers=10 xfers=5 timeout=60s; 
			print;
        		'$clock_ref';
    	    		set clock_out enable;
			trigger j51-1 rx '$trigger';
            		set biastee on;
            
			rx start;
			trigger '$triggerctrl' '$chain' '$fire';
			rx wait'
fi


# if Reception on Rx1 & Rx2 channels are required run the following
if [ "$tx_rx" = 'pass' ]
then 
	chain="rx";
	bladeRF-cli -d "*:serial=$sdr_serial" -v verbose -e '

			set frequency rx '$center_freq'M;
			set samplerate rx '$bw'M;
			set bandwidth rx '$bw'M;
			set agc off; 
			set gain rx1 '$rx1_gain'; 
			set gain rx2 '$rx2_gain'; 


			rx config file=/tmp/passive_'$test_id'.sc16q11 format=bin n='$cap_samps' channel=1,2 samples=134144 buffers=10 xfers=5 timeout=60s;

    			'$clock_ref';
    			set clock_out enable;
			trigger j51-1 rx '$trigger';
            		set biastee on;

            
			rx start;
            		trigger '$triggerctrl' '$chain' '$fire';
			rx wait'
fi
