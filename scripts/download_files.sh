#!/bin/bash
#$ -S /bin/bash
#$ -N downloaded_file
#$ -cwd
#$ -o /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -e /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -q short.q
#$ -l h_rt=24:00:00
#$ -pe thread 4
#$ -l h_vmem=2.75G

echo "JOB NAME: $JOB_NAME"
echo "JOB ID: $JOB_ID"
echo "QUEUE: $QUEUE"
echo "HOSTNAME: $HOSTNAME"
echo "SGE O WORKDIR: $SGE_O_WORKDIR"
echo "SGE TASK ID: $SGE_TASK_ID"
echo "NSLOTS: $NSLOTS"

echo Start at:`date`

conda activate pyega3_env
python3 /data1/scratch/pamesl/projet_cbf/scripts/selectFile.py /data1/scratch/pamesl/projet_cbf/Sample_File_SJCBF.map
conda deactivate

echo End at:`date`
