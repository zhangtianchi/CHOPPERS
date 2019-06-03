#!/bin/bash

rm log.*
#export I_MPI_COMPATIBILITY=4         #  Must

export MPI_NUM_THREADS=128




ulimit -m unlimited #memory
ulimit -s unlimited #stack
ulimit -v unlimited #virtual

for run in ccvt grid glass
do

input=/data/inspur_disk01/userdir/tczhang/work4/$run
output=/data/inspur_disk01/userdir/tczhang/work4/$run

fname=match-trackid_Subfind_Used


for snap in  ` seq 21 134`    # 14 31 57 84 134
#for snap in   134
do


mpirun -np $MPI_NUM_THREADS  ./HBT2LSubfind  $input $output  $snap $fname  >> log.$run

done
done


