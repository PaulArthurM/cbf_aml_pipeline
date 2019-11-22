#!/bin/bash
#$ -S /bin/bash
#$ -N oxog_picard
#$ -cwd
#$ -o /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -e /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -q short.q
#$ -l h_rt=24:00:00
#$ -pe thread 5
#$ -l h_vmem=2.75G

echo "JOB NAME: $JOB_NAME"
echo "JOB ID: $JOB_ID"
echo "QUEUE: $QUEUE"
echo "HOSTNAME: $HOSTNAME"
echo "SGE O WORKDIR: $SGE_O_WORKDIR"
echo "SGE TASK ID: $SGE_TASK_ID"
echo "NSLOTS: $NSLOTS"

echo Start at:`date`

conda activate RNAseq_pipeline_Gaetano

picard CollectOxoGMetrics \
      I=/data1/scratch/pamesl/projet_cbf/data/bam/SJCBF040_D.6.7_marked_duplicates_BQSR_merge.bam \
      O=/data1/scratch/pamesl/projet_cbf/data/oxog/oxoG_SJCBF040_metrics.txt \
      R=/data1/scratch/pamesl/projet_cbf/data/b37_data/GRCh37-lite.fa

conda deactivate

echo End at:`date`
