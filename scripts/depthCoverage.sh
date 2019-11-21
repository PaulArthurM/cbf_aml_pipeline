#!/bin/bash
#$ -S /bin/bash
#$ -N DepthOfCoverage_SJCBF040
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

conda activate gatk_env

gatk3 -T DepthOfCoverage \
  -omitBaseOutput \
  -omitLocusTable\
  -R /data1/scratch/pamesl/projet_cbf/data/b37_data/GRCh37-lite.fa \
  -I /data1/scratch/pamesl/projet_cbf/data/bam/SJCBF040_D.6.7_marked_duplicates_BQSR_merge.bam \
  -L /data1/scratch/pamesl/projet_cbf/data/b37_data/Broad.human.exome.b37.interval_list \
  -o SJCBF040_D.coverage \
  -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50

conda deactivate

echo End at:`date`
