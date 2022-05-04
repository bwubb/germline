


with open(config.get('project',{}).get('sample_list','sample.list'),'r') as s:
    SAMPLES=s.read().splitlines()

with open(config.get('project',{}).get('bam_table','bam.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

rule all:
    input:
        expand("data/work/{sample}/{library}/vardict.vcf.gz",sample=SAMPLES,library=config['resources']['targets_key'])

rule run_unpaired_vardict:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        "data/work/{sample}/{library}/vardict.vcf.gz"
    params:
        ref=config['reference']['fasta'],
        bed=config['resources']['targets_bed'],
        path="$HOME/software/VarDictJava/VarDict",
        AF_THR=0.01
    threads:
        4
    shell:
        "$HOME/software/VarDictJava/build/install/VarDict/bin/VarDict -th {threads} -G {params.ref} -f {params.AF_THR} -N {wildcards.sample} -b {input.bam} -z -c 1 -S 2 -E 3 -g 4 {params.bed} | {params.path}/teststrandbias.R | {params.path}/var2vcf_valid.pl -N {wildcards.sample} -f {params.AF_THR} | bgzip -c > {output}"
