#!/bin/sh

cap_samps=$1
pulses=$2
tx_gain=$3
rx_gain=$4
center_freq=$5
bw=$6
test_id=$7

echo "Running N0 bladeRF automated script"

#check if device is connected
bladeRF-cli -p 1>/dev/null 2>&1

if [ $? -ne 0 ]; then 
	echo "No bladeRF devices connected." >&2
	exit 1
fi


#Fetch device info

DEVICE_INFO=$(bladeRF-cli -e 'info' 2>/dev/null)

if [ -z "$DEVICE_INFO" ]; then 
	echo "Failed to open device and query info." >&2
	exit 1
fi 

#Read Device Serial Number 
DEVICE_SERIAL=$(echo "$DEVICE_INFO" | grep 'Serial #:' | sed -e 's/.*Serial #:\s\+//') 

#Read if device has FPGA loaded 
DEVICE_FPGA_LOADED=$(echo "$DEVICE_INFO" | grep "FPGA loaded:" | sed -e 's/.*FPGA loaded:\s\+//')

#Read device FPGA size
DEVICE_FPGA_SIZE=$(echo "$DEVICE_INFO" | grep "FPGA size:" | sed -e 's/.*FPGA size:\s\+//')

echo ""
echo "Found bladeRF"
echo "Serial #: $DEVICE_SERIAL"
echo "FPGA type: $DEVICE_FPGA_SIZE"
echo "FPGA Loaded: $DEVICE_FPGA_LOADED"
echo ""

now="$(date)"
echo "$now"
rxfilename=$now

echo "$rxfilename"


bladeRF-cli -d "*:serial=4eedd4abca9b48a0a49d2a469f09b719" -e 'set frequency rx '$center_freq'M;
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
