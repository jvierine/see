#!/bin/bash 

#sudo ntpdate ntp.uit.no

# AM 
freq1=5.5e6
rate=1e6

datestr=`date +%Y.%m.%d`
echo $datestr

dir1=/data0/see/${datestr}/see${rate}_${freq1}

mkdir -p ${dir1}

echo ${dir1}


nohup thor.py -m 192.168.10.2 -d "A:A" -c cha -f $freq1 -r $rate ${dir1}  > log.2 2>&1 &
