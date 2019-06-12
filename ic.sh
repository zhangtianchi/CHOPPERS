#!/bin/sh

rm nohup.*
rm -r Task*
#export I_MPI_COMPATIBILITY=4   

#if [ $# != 2  ] ; then 
#    echo "USAGE:    sh  $0  grid(glass, ccvt) beginsnap " 
#exit 1; 
#fi 

export MPI_NUM_THREADS=128

ulimit -m unlimited #memory
ulimit -s unlimited #stack
ulimit -v unlimited #virtual

files=128
haloname=orgin_halo
haloname1=halo

rm  /data/inspur_disk01/userdir/tczhang/work4/halo.*

for run in  grid #glass ccvt
#for run in   ccvt0 ccvt1 grid0 grid1 glass0 glass1
do

path=/data/inspur_disk01/userdir/tczhang/work4/$run

export RUN_NUM=$run
date >>nohup.$RUN_NUM

#for snap in  `seq 14 134`
for snap in  20 #31 57
#for snap in  134
do

fsnap=$(printf "%03d" $snap)
rm -r /data/inspur_disk01/userdir/tczhang/work4/$run/halo_$fsnap
rm -r /data/inspur_disk01/userdir/tczhang/work4/$run/orgin_halo_$fsnap

start=$(date +%s)
mpirun -np $MPI_NUM_THREADS   ./CHOPPERS  $path  $snap  $files $haloname    #1>>/data/inspur_disk01/userdir/tczhang/work4/halo_orgin.$RUN_NUM 2>&1 
end=$(date +%s)

#start1=$(date +%s)
#mpirun -np $MPI_NUM_THREADS   ./CHOPPERS_unbinding  $path  $snap  $files $haloname1    1>>/data/inspur_disk01/userdir/tczhang/work4/halo.$RUN_NUM 2>&1 
#end1=$(date +%s)

time=$(( $end - $start  ))
#time1=$(( $end1 - $start1  ))

echo "snap=" $snap ", no unbinding spend" $time " second, unbinding spend" $time1 " second!"  >>nohup.$RUN_NUM


done
done

#echo “ ALL DenFind Finished!!! ” | mail -s Denfind_ALL elsenhower@163.com












