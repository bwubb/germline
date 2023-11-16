


with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

def map_Workflow(wildcards):
    if 'targets_bedgz' in config['resources']:
        return rules.write_manta_wes.output
    return rules.write_manta_wgs.output

rule all_manta:
    input:
        expand("data/work/{lib}/{sample}/manta-wes/results/variants/candidateSV.vcf.gz",sample=SAMPLES,lib=config['resources']['targets_key'])

#Need better wes vs wgs
rule write_manta_wgs:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        temp("{work_dir}/{sample}/manta-wgs/runWorkflow.py")
    params:
        runDir="{work_dir}/{sample}/manta-wgs",
        reference=config['reference']['fasta']
    shell:
        "$HOME/software/manta/bin/configManta.py --bam {input} --referenceFasta {params.reference} --runDir {params.runDir}"

rule write_manta_wes:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        temp("{work_dir}/{sample}/manta-wes/runWorkflow.py")
    params:
        runDir="{work_dir}/{sample}/manta-wes",
        reference=config['reference']['fasta'],
        bedgz=config['resources'].get('targets_bedgz','')
    shell:
        "$HOME/software/manta/bin/configManta.py --bam {input.bam} --referenceFasta {params.reference} --callRegions {params.bedgz} --exome --runDir {params.runDir}"

#Looking for way to merge wes wgs runs.
rule run_manta:
    input:
        "{work_dir}/{sample}/manta-wes/runWorkflow.py"
        #map_Workflow
    output:
        "{work_dir}/{sample}/manta-wes/results/variants/candidateSmallIndels.vcf.gz",
        "{work_dir}/{sample}/manta-wes/results/variants/candidateSV.vcf.gz",
        "{work_dir}/{sample}/manta-wes/results/variants/diploidSV.vcf.gz",
    threads:
        4
    shell:
        "{input} -m local -j {threads}"
