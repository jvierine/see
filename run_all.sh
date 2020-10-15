#!/bin/bash 

sudo ntpdate ntp.uit.no

freq=4.04e6
rate=1e6

datestr=`date +%Y.%m.%d`
echo $datestr

dir0=/data0/${datestr}/test${rate}_${freq}
dir1=/data1/${datestr}/test${rate}_${freq}

mkdir -p ${dir0}
mkdir -p ${dir1}

echo ${dir0}
echo ${dir1}

nohup thor.py -m 192.168.10.3 -d "A:A A:B" -c cha,chb -f $freq -r $rate ${dir0}  > log.3 2>&1 &
nohup thor.py -m 192.168.10.4 -d "A:A A:B" -c chc,chd -f $freq -r $rate ${dir0}  > log.4 2>&1 &
nohup thor.py -m 192.168.10.5 -d "A:A A:B" -c che,chf -f $freq -r $rate ${dir1}  > log.5 2>&1 &
nohup thor.py -m 192.168.10.6 -d "A:A A:B" -c chg,chh -f $freq -r $rate ${dir1}  > log.6 2>&1 &

