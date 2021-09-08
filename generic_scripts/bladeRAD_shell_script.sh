#!/bin/sh

test_id=$1
cap_samps=$2
pulses=$3
tx_gain=$4
rx_gain=$5
center_freq=$6
bw=$7
sdr_serial=$8
trigger=$9
clock=${10}
tx_rx=${11}

#Create trigger message
if [ $trigger='master' ]; then triggerctrl="trigger j51-1 $tx_rx fire";
else triggerctrl=''
fi	

echo $triggerctrl

#creatr clock ref message
if [ $clock=1 ]; then clock_ref='set clock_sel external';
elif [ $clock=2 ]; then clock_ref='set clock_out enable';
else clock_ref='set clock_ref enable'; fi


# if Transmission required run the following
if [ $tx_rx='tx' ]
then 
	bladeRF-cli -d "*:serial=$sdr_serial" -e ' set frequency tx '$center_freq'M;
			set samplerate tx '$bw'M;
			set bandwidth tx '$bw'M;
			set gain tx1 '$tx_gain';

			tx config file=/tmp/chirp.sc16q11 format=bin samples=16384 buffers=32 repeat='$pulses' timeout=30s;
	
			trigger j51-1 tx '$trigger';

			tx start;
			'$triggerctrl';
			tx wait;
			print'
fi


