#/bin/sh
#BSUB -J ciCor[1-100]
#BSUB -o logs/ciCor.out%J
#BSUB -e logs/ciCor.err%J
#BSUB -M 1000
#BSUB -R "select[mem>1000] rusage[mem=1000]" 
R --no-save < ciCor.R