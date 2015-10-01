#/bin/sh
#BSUB -J imputation[1-1540]
#BSUB -o logs/imputation.out%J
#BSUB -e logs/imputation.err%J
#BSUB -M 1000
#BSUB -R "select[mem>1000] rusage[mem=1000]" 
R --no-save < imputation.R