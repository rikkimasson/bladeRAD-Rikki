#!/bin/sh

test_id='test1'
cap_samps='400M'
bw=60000000



bladeRF-cli -v critical -e '
			
			set samplerate rx '$bw';
			set gpio 0x257;
			
			rx config file=/tmp/sample_rate_test.sc16q11 format=bin n='$cap_samps' samples=134144 buffers=32 xfers=16 timeout=30s; 
            
            set frequency tx 70M;
			set samplerate tx '$bw';
			set gain tx1 0;

			tx config file=/tmp/tone format=bin samples=16384 buffers=32 repeat=292968 timeout=60s;
	
            
			trigger j51-1 rx master;
            trigger j51-1 tx slave;

            
            
			rx start;
            tx start;
			trigger j51-1 rx fire;
			rx wait;
            tx wait'
			
			