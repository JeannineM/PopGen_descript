#!/bin/bash
#SBATCH -J simul_model1        	## name
#SBATCH -N 1                 	## 1 node
#SBATCH --ntasks=16          	## number of tasks per node  
#SBATCH --mail-type=END
#SBATCH --mail-user=blanckaert.a@gmail.com


mySTART=`pwd`;

cd $SCRATCH
rm -rf alextmp
mkdir alextmp
cd alextmp

module purge
module load gcc/4.8.2
module load intel-mpi/5
module load intel/14.0.2
module load gsl/1.16
module load python/2.7
module load R/3.1.1 			## maybe not the proper syntax, I will check later


cp $mySTART/simulator_cluster.R ./
cp $mySTART/make_bash.py ./
cp $mySTART/run_FD_4pop_cluster.R ./
cp $mySTART/jsfs_identity.csv ./

python make_bash.py 		## this create 16 different bash files (name:bash_script_n .sh)

chmod +x bash_script_n*

for ((i=1;i<=16;i++)); do ./bash_script_n${i}.sh & done ## run all 16 jobs, one per core , each job does a certain number of simulations
wait

cp * $mySTART/result	## copy back the results
cd ..
rm -rf alextmp
