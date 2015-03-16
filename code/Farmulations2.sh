#/bin/sh
#BSUB -J Farmulations[1-1000]
#BSUB -o logs/Farmulations.out%J
#BSUB -e logs/Farmulations.err%J
#BSUB -M 1000
#BSUB -R "select[mem>1000] rusage[mem=1000]" 
export NJOBS=100
R --no-save < Farmulations2.R
