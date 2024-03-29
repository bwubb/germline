
import datetime
import sys
import os
import math
import errno

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno==errno.EEXIST and os.path.isdir(path):
            pass


def bam_input(wildcards):
    return BAMS[wildcards.sample]

def gvcf_input(wildcards):
    return ' '.join([f'-V data/work/{wildcards.reference}/{sample}/gatk/haplotype.g.vcf.gz' for sample in SAMPLES])

def cohort_gvcf_input(wildcards):
    return ' '.join([f"-V data/work/{wildcards.reference}/{wildcards.project}/gatk/cohort{n}.haplotype.g.vcf.gz" for n in range(1,cohort_num+1)])

def chr_gvcf(wildcards):
    return expand("data/work/{reference}/{sample}/gatk/chr{chr}.g.vcf.gz",sample=SAMPLES,chr=wildcards.chr,reference=wildcards.reference)

def chr_cohort(wildcards):
    return [f"data/work/{wildcards.reference}/{wildcards.project}/gatk/cohort{n}_chr{wildcards.chr}.g.vcf.gz" for n in range(1,cohort_num+1)]

def modrange(first, last, modulus, step = 1):
    for i in range(first, last, step):
        yield i % modulus

with open(config['project']['sample_list'],'r') as infile:
    SAMPLES=infile.read().splitlines()
    for sample in SAMPLES:
        mkdir_p('logs/cluster/%s' % sample)

with open(config['project']['bam_list'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

mkdir_p(f'logs/cluster/{config["project"]["name"]}')

#***100 or config
cohort_num=math.ceil(len(SAMPLES)/75)
CHR=list(range(1,23))+['X']

wildcard_constraints:
    chr="[0-9]{1,2}|[X-Y]{1}",
    #work_dir=f"data/work/{config['resources']['targets_key']}",
    #results_dir=f"data/final/{config['project']['name']}"


rule all:
    input: expand("data/work/{reference}/{project}/gatk/haplotype.chr{chr}.vcf.gz",reference="GRCh37",project=config['project']['name'],chr=CHR)

rule all_gvcf:
    input:
        expand('data/work/{reference}/{sample}/gatk/haplotype.g.vcf.gz',reference=config['reference']['key'],sample=SAMPLES)

rule HaplotypeCaller_GVCF:
    input:
        bam_input
    output:
        "data/work/{reference}/{sample}/gatk/haplotype.g.vcf.gz"
    params:
        reference=config['reference']['fasta'],
        memory='10240m'
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T HaplotypeCaller -I {input} --emitRefConfidence GVCF -o {output} -G Standard -G AS_Standard"

rule HaplotypeCaller_byChr:
    input:
        bam_input
    output:
        "data/work/{reference}/{sample}/gatk/chr{chr}.g.vcf.gz"
    params:
        reference=config['reference']['fasta'],
        chr=lambda wildcards: wildcards.chr
    shell:
        "java -Xmx10g -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T HaplotypeCaller -I {input} -L {params.chr} --emitRefConfidence GVCF -o {output} -G Standard -G AS_Standard"

rule Make_Cohorts_byChr:
    input:
        chr_gvcf
    output:
        expand("data/work/{reference}/{project}/gatk/chr{chr}.cohort{n}.list",reference=config['reference']['key'],project=config['project']['name'],chr=CHR,n=list(range(1,cohort_num+1)))
    script:
        'make_cohorts.py'

rule CombineGVCFs_Cohorts_byChr:
    input:
        "data/work/{reference}/{project}/gatk/chr{chr}.cohort{n}.list"
    output:
        "data/work/{reference}/{project}/gatk/chr{chr}.cohort{n}.g.vcf.gz"
    params:
        reference=config['reference']['fasta'],
        #variants=gvcf_input,
        memory='32g'
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T CombineGVCFs -V {input} -o {output}"

rule GenotypeGVCFs_byChr:
    input:
        chr_cohort
    output:
        "data/work/{reference}/{project}/gatk/haplotype.chr{chr}.vcf.gz"
    params:
        reference=config['reference']['fasta'],
        memory='32g'
    run:
        #Could also use unsorted input
        GVCFS=[f'-V {i}' for i in input]
        shell(f"java -Xmx{params.memory} -jar $HOME/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -R {params.reference} -T GenotypeGVCFs {' '.join(GVCFS)} --disable_auto_index_creation_and_locking_when_reading_rods -o {output}")
