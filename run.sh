#!/bin/sh


export MPI_NUM_THREADS=16

ulimit -m unlimited #memory
ulimit -s unlimited #stack
ulimit -v unlimited #virtual

files=128

haloname=halo

path=/home/tczhang/work1

snap=100


mpirun -np $MPI_NUM_THREADS   ./CHOPPERS  $path  $snap  $files $haloname    1>>$path/log.txt 2>&1 













