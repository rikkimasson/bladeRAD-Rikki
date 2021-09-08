#!/bin/sh

cap_samps=$1
pulses=$2
tx_gain=$3
rx_gain=$4
center_freq=$5
bw=$6
test_id=$7

sdr_serial=$8
trigger=$9




bladeRF-cli -d "*:serial='$sdr_serial'" -e 'set frequency rx '$center_freq'M;
		set samplerate rx '$bw'M;
		set bandwidth rx '$bw'M;
		set agc off; 
		rx config file=/tmp/'$test_id'.sc16q11 format=bin n='$cap_samps' samples=16384 buffers=32 xfers=16 timeout=30s;
		set gain rx1'$rx_gain' ; 

		set frequency tx '$center_freq'M;
		set samplerate tx '$bw'M;
		set bandwidth tx '$bw'M;
		set gain tx1 '$tx_gain';
		set gain tx2 -20;
		set gain tx2 dsa -89;

		tx config file=/tmp/40MHz_chirp.sc16q11 format=bin samples=65536 buffers=32 repeat='$pulses' timeout=30s;

		
		trigger j51-1 tx master;
		trigger j51-1 rx slave;
		
		print;

	
		tx start;
		rx start;

		trigger j51-1 tx fire;
	 
		tx wait;
		rx wait'

 mv /tmp/$test_id.sc16q11 /home/piers/Documents/bladeRF_Code/FMCW/Experiments
