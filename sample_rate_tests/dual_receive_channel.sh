#!/bin/sh

test_id='test1'
cap_samps='400M'
bw=20000000



bladeRF-cli -e '
			
			set samplerate rx '$bw';
			set gpio 0x257;
			
			rx config file=/tmp/sample_rate_test.sc16q11 channel=1,2 format=bin n='$cap_samps' samples=134144 buffers=32 xfers=16 timeout=30s; 
		
			trigger j51-1 rx master;
            
            
			rx start;
			trigger j51-1 rx fire;
			rx wait'
			
			