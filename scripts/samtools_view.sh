#!/bin/bash
#$ -S /bin/bash
#$ -N samtools_view
#$ -cwd
#$ -o /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -e /data1/scratch/pamesl/projet_cbf/stdoe_sge
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

echo Start at:`date`

bam_file=/data1/scratch/pamesl/projet_cbf/data/bam/samples/$1

echo Work on ${bam_file} file.

conda activate samtools_env
samtools view -F 0x400 $bam_file
conda deactivate

echo End at:`date`
