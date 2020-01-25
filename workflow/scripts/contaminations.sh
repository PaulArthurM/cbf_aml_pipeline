#!/bin/bash
#$ -S /bin/bash
#$ -N contaminationVCF
#$ -cwd
#$ -o /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -e /data1/scratch/pamesl/projet_cbf/stdoe_sge
#$ -q short.q
#$ -l h_rt=01:00:00
#$ -pe thread 8
#$ -l h_vmem=2.75G

echo "JOB NAME: $JOB_NAME"
echo "JOB ID: $JOB_ID"
echo "QUEUE: $QUEUE"
echo "HOSTNAME: $HOSTNAME"
echo "SGE O WORKDIR: $SGE_O_WORKDIR"
echo "SGE TASK ID: $SGE_TASK_ID"
echo "NSLOTS: $NSLOTS"

echo Start at:`date`

conda activate gatk4_4.1.2.0_env
gatk SelectVariants -V /data1/scratch/pamesl/projet_cbf/data/GnomAD/af-only-gnomad.raw.sites.b37.vcf --select 'AF > 0.05' -O /data1/scratch/pamesl/projet_cbf/data/vcf/variants_for_contamination.vcf
conda deactivate

echo End at:`date`
