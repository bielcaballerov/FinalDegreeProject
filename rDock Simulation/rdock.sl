#!/bin/bash
#SBATCH -J rdock
#SBATCH --output=rdock.out
#SBATCH --error=rdock.err
#SBATCH --ntasks=16
#SBATCH --time=00-02:00:00

module purge
module load intel mkl impi gcc
module load impi
module load boost/1.64.0
module load rdock/2013.1

#Create the cavity where the ligands are going to be bound
rbcavity -r standard_rDock.prm -was

#Docking simulation
rbdock -r standard_rDock.prm -p dock.prm -n 5 -i ligands.sdf -o ligands_out

#Select the best poses from previous step, for each ligand
sdsort -n -s -fSCORE ligands_out.sd | sdfilter -f'$_COUNT == 1' > bestposes.sd

#Make a report to obtain the rDock Score
sdreport -t bestposes.sd | awk '{print $2,$3,$4,$5,$6,$7}' > bestposes.csv
