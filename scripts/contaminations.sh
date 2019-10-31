#!/bin/bash
#$ -S /bin/bash
#$ -N samtools_index
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

conda activate gatk4_4.1.2.0_env
gatk SelectVariants -V data/GnomAD/af-only-gnomad.raw.sites.b37.vcf -L data/b37_data/Broad.human.exome.b37.interval_list --select "AF > 0.05" -O data/vcf/^Criants_for_contamination.vcf
conda deactivate

echo End at:`date`
