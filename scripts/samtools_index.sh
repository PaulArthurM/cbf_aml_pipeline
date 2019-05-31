#!/bin/bash
#$ -S /bin/bash
#$ -N samtools_index
#$ -cwd
#$ -o /data1/scratch/pamesl/stdoe_sge
#$ -e /data1/scratch/pamesl/stdoe_sge
#$ -q short.q
#$ -l h_rt=01:00:00
#$ -pe thread 1
#$ -l h_vmem=2.75G

echo "JOB NAME: $JOB_NAME"
echo "JOB ID: $JOB_ID"
echo "QUEUE: $QUEUE"
echo "HOSTNAME: $HOSTNAME"
echo "SGE O WORKDIR: $SGE_O_WORKDIR"
echo "SGE TASK ID: $SGE_TASK_ID"
echo "NSLOTS: $NSLOTS"

IN_BAM=/data1/scratch/pamesl/projet_cbf/data/bam/${1}

conda activate samtools_env
samtools index -b ${IN_BAM}.bam ${IN_BAM}.bai 
conda deactivate samtools_env
