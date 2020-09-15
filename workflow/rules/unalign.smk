
input:
    "results/preprocessing/{sample}_{type}.bam"
output:
    "results/preprocessing/{sample}_{type}.unaligned.bam"
log:
    "logs/revertsam/{sample}.log"
conda:
    "../envs/picard.yaml"
shell:
    "picard RevertSam \
        I={input} \
        O={output} \
        SANITIZE=true \
        MAX_DISCARD_FRACTION=0.005 \
        ATTRIBUTE_TO_CLEAR=XT \
        ATTRIBUTE_TO_CLEAR=XN \
        ATTRIBUTE_TO_CLEAR=AS \
        ATTRIBUTE_TO_CLEAR=OC \
        ATTRIBUTE_TO_CLEAR=OP \
        SORT_ORDER=queryname \
        RESTORE_ORIGINAL_QUALITIES=true \
        REMOVE_DUPLICATE_INFORMATION=true \
        REMOVE_ALIGNMENT_INFORMATION=true"
