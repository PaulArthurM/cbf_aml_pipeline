#!/bin/bash
#$ -S /bin/bash
#$ -N samtools_depth
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

echo Work on: $1

module load samtools/1.2
samtools depth -a $1 | awk '{c++;s+=$3}END{print s/c}'
module unload samtools/1.2

echo End at:`date`
