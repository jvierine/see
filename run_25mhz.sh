#!/usr/bin/bash

thor.py -m 192.168.10.3 -d A:A -c cha -f 12.5e6 -r 25e6 /dev/shm/hf25 &
sleep 5
while true; do
    rsync -av --remove-source-files --exclude=tmp* --progress /dev/shm/hf25/cha /data_out/hf25/
    sleep 1
done
