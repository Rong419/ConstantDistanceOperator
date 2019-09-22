#!/bin/bash -e
#SBATCH -J SummaryTree
#SBATCH -A nesi00390
#SBATCH --time=00:50:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --hint=nomultithread
#SBATCH -D /nesi/project/nesi00390/rong/operator/calibrated/summarytree
#SBATCH -o s_Output_FILE.txt
#SBATCH -e s_Error_FILE.txt

module load beagle-lib/3.0.1-gimkl-2017a
module load Java/1.8.0_144

srun java -Xmx4096m -Djava.library.path=$BEAGLE_LIB_PATH -jar /nesi/project/nesi00390/rong/operator/TreeAnnotator.jar -heights mean -b 15 /nesi/project/nesi00390/rong/operator/calibrated/output/FILE.trees /nesi/project/nesi00390/rong/operator/calibrated/summarytree/s_FILE.trees

