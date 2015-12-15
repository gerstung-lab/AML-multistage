#/bin/sh
#BSUB -J leaveOneOut[1-1540]
#BSUB -o logs/leaveOneOut.out%J
#BSUB -e logs/leaveOneOut.err%J
#BSUB -M 1000
#BSUB -R "select[mem>1000] rusage[mem=1000]" 
R --no-save < leaveOneOut.R