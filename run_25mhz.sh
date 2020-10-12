#!/usr/bin/bash
#
# Start a ram disk ringbuffer with linux and digital rf
# This is for realtime processing
#
drf ringbuffer -z 4000MB /dev/shm/hf25 -p 2 &
thor.py -m 192.168.10.3 -d A:A -c cha -f 12.5e6 -r 25e6 /dev/shm/hf25 
