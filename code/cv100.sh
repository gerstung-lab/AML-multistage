#/bin/sh
#BSUB -J cv100[1-100]
#BSUB -o logs/cv100.out%J
#BSUB -e logs/cv100.err%J
#BSUB -M 1000
#BSUB -R "select[mem>1000] rusage[mem=1000]" 
R --no-save < cv100.R