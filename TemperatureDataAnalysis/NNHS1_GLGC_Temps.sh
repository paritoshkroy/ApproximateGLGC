#!/bin/bash
#SBATCH --account=def-aschmidt  # replace this with your own account
#SBATCH --ntasks=16              # number of processes
#SBATCH --mem-per-cpu=16000M      # memory; default unit is megabytes
#SBATCH --time=48:00:00         # time (HH:MM:SS)
#SBATCH --output=/home/pkroy/projects/def-aschmidt/pkroy/ApproximateGLGC/TemperatureDataAnalysis/%x-%j.out

# Modules
module load gcc/11.3.0
module load r/4.2.1



# Export the nodes names. 
# If all processes are allocated on the same node, NODESLIST contains : node1 node1 node1 node1
# Cut the domain name and keep only the node name
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))
R -f NNHS1_GLGC_Temps.R
