#!/bin/bash
#$ -S /bin/bash
#$ -N unzip
#$ -cwd
#$ -o /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -e /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -q short.q
#$ -l h_rt=01:00:00
#$ -pe thread 4
#$ -l h_vmem=11G

echo "JOB NAME: $JOB_NAME"
echo "JOB ID: $JOB_ID"
echo "QUEUE: $QUEUE"
echo "HOSTNAME: $HOSTNAME"
echo "SGE O WORKDIR: $SGE_O_WORKDIR"
echo "SGE TASK ID: $SGE_TASK_ID"
echo "NSLOTS: $NSLOTS"

echo Start at:`date`

gzip -d ${1}

echo End at:`date`
