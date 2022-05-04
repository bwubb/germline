#! snakemake
import os


with open(config.get('project',{}).get('sample_list','samples.list'),'r') as infile:#sample_list in config, default is samples.list
    SAMPLES=infile.read().splitlines()

os.makedirs(f'logs/cluster/{config["project"]["name"]}',exist_ok=True)
for sample in SAMPLES:
    os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())
    #I should really check to make sure they exist...

chr_list=list(range(1,23))+['X','Y']
chr_dict={}
if config['reference']['key'] in ['GRCh37','GRCh38']:
    for c in chr_list:
        chr_dict[f'chr{c}']=f'{c}'
elif config['reference']['key'] in ['hg19','hg38']:
    for c in chr_list:
        chr_dict[f'chr{c}']=f'chr{c}'
else:
    print('Unknown reference key. Please include ref.dict')
    raise

def bam_input(wildcards):
    return BAMS[wildcards.sample]

def chr_gvcf(wildcards):
    return expand("data/work/{reference}/{sample}/gatk/{chr}.g.vcf.gz",sample=SAMPLES,chr=wildcards.chr,reference=wildcards.reference)

def gvcf_input(wildcards):
    return ' '.join([f"-V data/work/{wildcards.reference}/{sample}/gatk/{wildcards.chr}.g.vcf.gz" for sample in SAMPLES])

def chr_func(wildcards):
    return chr_dict[wildcards.chr]

#localrules:

rule byChr:
    input:
        expand("data/work/{reference}/{project}/gatk/haplotype.{chr}.vcf.gz",reference=config['reference']['key'],project=config['project']['name'],chr=chr_dict.keys())

rule GATK4_HaplotypeCaller_byChr:
    input:
        bam_input
    output:
        "data/work/{reference}/{sample}/gatk/{chr}.g.vcf.gz"
    params:
        reference=config['reference']['fasta'],
        chr=chr_func,
        memory='10g'
    shell:
        "gatk --java-options '-Xmx{params.memory}' HaplotypeCaller -I {input} -O {output} -R {params.reference} --emit-ref-confidence GVCF -L {params.chr}"

rule GATK4_CombineGVCFs_byChr:
    input:
        chr_gvcf
    output:
        "data/work/{reference}/{project}/gatk/{chr}.g.vcf.gz"
    params:
        reference=config['reference']['fasta'],
        variants=gvcf_input,
        memory='30g'
    shell:
        "gatk --java-options '-Xmx{params.memory}' CombineGVCFs {params.variants} -O {output} -R {params.reference}"

rule GATK4_GenotypeGVCFs_byChr:
    input:
        "data/work/{reference}/{project}/gatk/{chr}.g.vcf.gz"
    output:
        "data/work/{reference}/{project}/gatk/haplotype.{chr}.vcf.gz"
    params:
        reference=config['reference']['fasta'],
        memory="30g"
    shell:
        "gatk --java-options '-Xmx{params.memory}' GenotypeGVCFs -R {params.reference} -V {input} -O {output}"

#No clue what are currect germline filter recommendations.
