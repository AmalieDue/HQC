#BSUB -J SageMath
#BSUB -q hpc
#BSUB -n 4
#BSUB -R "rusage[mem=2GB]"
#BSUB -R "span[hosts=1]"
#BSUB -R "select[avx2]"
#BSUB -W 3:00

/appl/SageMath/9.4/sage simulation.py
