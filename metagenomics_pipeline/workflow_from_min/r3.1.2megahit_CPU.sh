#!/bin/bash
pid=334989  #获取进程pid
echo $pid
interval=10  #设置采集间隔
while true
do
    echo $(date +"%y-%m-%d %H:%M:%S") >> test.txt
    cpu=`qstat -j $pid |grep 'usage'`    #获取cpu占用
    echo "Cpu:  $cpu" >> test.txt
    echo $blank >> test.txt
    sleep $interval
done
