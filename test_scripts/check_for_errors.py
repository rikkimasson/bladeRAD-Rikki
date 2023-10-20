#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 13:48:23 2022

@author: piers
"""

#!env python3
import sys
import array

if len(sys.argv) != 2:
    print('Usage: ' + sys.argv[0] + ': <data file>')
    sys.exit(1)

with open(sys.argv[1], 'rb') as in_file:
    data = in_file.read()
    count = array.array('I')
    count.fromstring(data)
    
    fail = 0
    curr = count[0]
    for i in range(1, len(count)):
        exp = curr + 1
    
        if count[i] != exp:
            print('[' + str(i) + '] = ' + str(count[i]) + ', Expected ' + str(exp) +
                  ', Gap = ' + str(exp - count[i]))
            fail += 1
    
        curr = count[i]
    
    print('Number of gaps:' + str(fail))