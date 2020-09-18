#!/bin/bash 

sudo ntpdate ntp.uit.no

freq0=4.15e6
freq1=4.15e6
freq2=6.225e6
freq3=8.3e6
rate=1e6

datestr=`date +%Y.%m.%d`
echo $datestr

dir0=/data0/${datestr}/test${rate}_${freq0}
dir1=/data0/${datestr}/test${rate}_${freq1}
dir2=/data1/${datestr}/test${rate}_${freq2}
dir3=/data1/${datestr}/test${rate}_${freq3}


mkdir -p ${dir0}
mkdir -p ${dir1}
mkdir -p ${dir2}
mkdir -p ${dir3}

echo ${dir0}
echo ${dir1}
echo ${dir2}
echo ${dir3}

nohup thor.py -m 192.168.10.3 -d "A:A A:B" -c cha,chb -f $freq0 -r $rate ${dir0}  > log.3 2>&1 &
nohup thor.py -m 192.168.10.4 -d "A:A A:B" -c chc,chd -f $freq1 -r $rate ${dir1}  > log.4 2>&1 &
nohup thor.py -m 192.168.10.5 -d "A:A A:B" -c che,chf -f $freq2 -r $rate ${dir2}  > log.5 2>&1 &
nohup thor.py -m 192.168.10.6 -d "A:A A:B" -c chg,chh -f $freq3 -r $rate ${dir3}  > log.6 2>&1 &
