#!/bin/bash
#$ -S /bin/bash
#$ -N SelectVariants
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

echo Variant file: ${1}
echo Output VCF: ${2}

conda activate gatk4_4.1.2.0_env
gatk SelectVariants \
    -V ${1} \
    -select-type SNP -restrict-alleles-to BIALLELIC \
    -O ${2}.vcf.gz
conda deactivate

echo End at:`date`
