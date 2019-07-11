#!/bin/bash
#$ -S /bin/bash
#$ -N IndexFeatures
#$ -cwd
#$ -o /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -e /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -q short.q
#$ -l h_rt=24:00:00
#$ -pe thread 6
#$ -l h_vmem=16.5G

echo "JOB NAME: $JOB_NAME"
echo "JOB ID: $JOB_ID"
echo "QUEUE: $QUEUE"
echo "HOSTNAME: $HOSTNAME"
echo "SGE O WORKDIR: $SGE_O_WORKDIR"
echo "SGE TASK ID: $SGE_TASK_ID"
echo "NSLOTS: $NSLOTS"

echo Start at:`date`

conda activate gatk4_4.1.2.0_env
gatk IndexFeatureFile \
     -F ${1}
conda deactivate

echo End at:`date`
