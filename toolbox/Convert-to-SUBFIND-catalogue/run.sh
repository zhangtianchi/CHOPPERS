#!/bin/bash

rm log.*

export MPI_NUM_THREADS=128




ulimit -m unlimited #memory
ulimit -s unlimited #stack
ulimit -v unlimited #virtual

for run in c cc ccc
do

input=/data/inspur_disk01/userdir/tczhang/work4/$run
output=/data/inspur_disk01/userdir/tczhang/work4/$run

fname=match-trackid_Subfind_Used


for snap in  ` seq 21 134`    # 14 31 57 84 134
do


mpirun -np $MPI_NUM_THREADS  ./HBT2LSubfind  $input $output  $snap $fname  >> log.$run

done
done


