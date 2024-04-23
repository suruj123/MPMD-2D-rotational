#!/bin/bash
##declare a name for this job to be sample_job
#PBS -N Ckh_k2.0_p5_0.5_th_g450_rout450
##request the parallel queue for this job
#PBS -q mediumq
##request a total of 28 processors for this job (1 nodes and 24 processors per node)
#PBS -l select=15:ncpus=40
##Request to run for specified hrs:min:secs
####PBS -l walltime=48:00:00
##combine PBS standard output and error files
#PBS -j oe
#PBS -V
##mail is sent to you when the job starts and when it terminates or aborts
#PBS -m bea
##specify your email address
###PBS -M suruj.kalita@ipr.res.in
##change to the directory where you submitted the job
module load intel-2018
#module load gsl/2.3
cd $PBS_O_WORKDIR
PREFIX=Ckh_k2.0_p5_0.5_th_g450_rout450
#include the full path to the name of your MPI program
mpiexec.hydra -n 600 /scratch/scratch_run/suruj.kalita/MPMD_2D_brute_mpi_circular_kh_2/source/main $PREFIX
mkdir ../output/$PREFIX
mv ../output/$PREFIX.* ../output/$PREFIX
exit 0
