def get_sample_name(sample):
    return re.match("(.+?)\.bam$", sample).group(1)


def get_id(sample):
    return re.search("-(\w+)\.", sample).group(1)


def get_lane(sample):
    return re.search("\.(\d+)$", sample).group(1)


def getLane(sample):
    return re.search("\.(\d)\.bam$", sample).group(1)


def get_tumour_lanes(sample):
    return ".".join([get_lane(lane) for lane in sample['D']])


def get_normal_lanes(sample):
    return ".".join([get_lane(lane) for lane in sample['G']])




def get_input():
    wanted_input = []
    # Load json configuration file
    CONFIG_JSON = json.load(open(config["SAMPLES"]))
    SAMPLES = CONFIG_JSON['samples']
    wanted_input.expend(expand(config["PROJECT_DIR"] + "data/bam/{sample}_{type}_marked_duplicates_BQSR_merge.bam", sample=SAMPLES, type=['G', 'D']))
    wanted_input.expend(expand(config["PROJECT_DIR"] + "data/bam/{sample}_{type}_marked_duplicates_BQSR_merge.bai", sample=SAMPLES, type=['G', 'D']))
    wanted_input.expend(expand(config["PROJECT_DIR"] + "{sample}_{type}_marked_duplicates_BQSR_merge_fastqc.html", sample=SAMPLES, type=['G', 'D']))
    wanted_input.expend(expand(config["PROJECT_DIR"] + "data/vcf/{sample}_{type}_marked_duplicates_BQSR_merge_for_pon.vcf.gz", sample=SAMPLES, type=['G', 'D']))
    wanted_input.expend(expand(config["PROJECT_DIR"] + "results/variantCalling/strelka/{sample}_vcf", sample=SAMPLES))
    wanted_input.expend(expand(config["PROJECT_DIR"] + "data/vcf/{sample}_somatic.vcf.gz", sample=SAMPLES))
    wanted_input.expend(expand(config["PROJECT_DIR"] + "data/vcf/filtered/{sample}_somatic_filtered.vcf.gz", sample=SAMPLES))
    wanted_input.expend(expand(config["PROJECT_DIR"] + "data/vcf/annotated/{sample}.avinput", sample=SAMPLES))
    wanted_input.expend(expand(config["PROJECT_DIR"] + "data/vcf/filtered/{sample}_somatic_filtered_pass.vcf", sample=SAMPLES))
    wanted_input.expend(expand(config["PROJECT_DIR"] + "results/variantCalling/freebayes/{sample}.vcf", sample=SAMPLES))

    # for SAMPLE in SAMPLES:
    #     for TYPE in SAMPLES[SAMPLE]:
    #         LANES = SAMPLES[SAMPLE][TYPE]
            #file_1 = "{project_dir}data/bam/{sample}_{type}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane=get_lane(LANES[0]))
            #file_2 = "{project_dir}data/bam/{sample}_{type}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane=get_lane(LANES[1]))
            #if ( (len(LANES)==2) and (os.path.isfile(file_1)) and (os.path.isfile(file_2)) ):
            # wanted_input.append("{project_dir}data/bam/{sample}_{type}_marked_duplicates_BQSR_merge.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
            # wanted_input.append("{project_dir}data/bam/{sample}_{type}_marked_duplicates_BQSR_merge.bai".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
            # wanted_input.append("{project_dir}{fastq_dir}{sample}_{type}_marked_duplicates_BQSR_merge_fastqc.html".format(fastq_dir=config["FASTQC"]["DIR"], project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
            # wanted_input.append("{project_dir}data/vcf/{sample}_{type}_marked_duplicates_BQSR_merge_for_pon.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
            # wanted_input.append("results/variantCalling/strelka/{sample}_vcf")
            #VCF_IDX.append("{project_dir}data/vcf/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge_for_pon.vcf.idx".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
            # vcf_somatic = "{project_dir}data/vcf/{sample}_somatic.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=get_normal_lanes(SAMPLES[SAMPLE]), lanes_tumour=get_tumour_lanes(SAMPLES[SAMPLE]))
            # vcf_filtered = "{project_dir}data/vcf/filtered/{sample}_somatic_filtered.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=get_normal_lanes(SAMPLES[SAMPLE]), lanes_tumour=get_tumour_lanes(SAMPLES[SAMPLE]))
            # vcf_annotated = "{project_dir}data/vcf/annotated/{sample}.avinput".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=get_normal_lanes(SAMPLES[SAMPLE]), lanes_tumour=get_tumour_lanes(SAMPLES[SAMPLE]))
            # pass_vcf="{project_dir}data/vcf/filtered/{sample}_somatic_filtered_pass.vcf".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=get_normal_lanes(SAMPLES[SAMPLE]), lanes_tumour=get_tumour_lanes(SAMPLES[SAMPLE]))
            # if not vcf_somatic in VCF_SOMATIC:
            #     wanted_input.append(vcf_somatic)
            # if not vcf_filtered in VCF_FILERED:
            #     wanted_input.append(vcf_filtered)
            # if not vcf_annotated in VCF_ANNOTATED:
            #     wanted_input.append(vcf_annotated)
            #     wanted_input.append(pass_vcf)
            # elif (len(LANES)==3):
            #     file_3 = "{project_dir}data/bam/{sample}_{type}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane=get_lane(LANES[2]))
            #     if (os.path.isfile(file_3)):
            #         wanted_input.append("{project_dir}data/bam/{sample}_{type}_marked_duplicates_BQSR_merge.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
            #         wanted_input.append("{project_dir}data/bam/{sample}_{type}_marked_duplicates_BQSR_merge.bai".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
            #         wanted_input.append("{project_dir}{fastq_dir}{sample}_{type}_marked_duplicates_BQSR_merge_fastqc.html".format(fastq_dir=config["FASTQC"]["DIR"], project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
            #         wanted_input.append("{project_dir}data/vcf/{sample}_{type}_marked_duplicates_BQSR_merge_for_pon.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
            #         #VCF_IDX.append("{project_dir}data/vcf/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge_for_pon.vcf.idx".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
            #         vcf_somatic = "{project_dir}data/vcf/{sample}_somatic.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=get_normal_lanes(SAMPLES[SAMPLE]), lanes_tumour=get_tumour_lanes(SAMPLES[SAMPLE]))
            #         vcf_filtered = "{project_dir}data/vcf/filtered/{sample}_somatic_filtered.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=get_normal_lanes(SAMPLES[SAMPLE]), lanes_tumour=get_tumour_lanes(SAMPLES[SAMPLE]))
            #         vcf_annotated = "{project_dir}data/vcf/annotated/{sample}.avinput".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=get_normal_lanes(SAMPLES[SAMPLE]), lanes_tumour=get_tumour_lanes(SAMPLES[SAMPLE]))
            #         pass_vcf="{project_dir}data/vcf/filtered/{sample}_somatic_filtered_pass.vcf".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=get_normal_lanes(SAMPLES[SAMPLE]), lanes_tumour=get_tumour_lanes(SAMPLES[SAMPLE]))
            #         if not vcf_somatic in VCF_SOMATIC:
            #             wanted_input.append(vcf_somatic)
            #         if not vcf_filtered in VCF_FILERED:
            #             wanted_input.append(vcf_filtered)
            #         if not vcf_annotated in VCF_ANNOTATED:
            #             wanted_input.append(vcf_annotated)
            #             wanted_input.append(pass_vcf)
    return wanted_input
