


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

rule brca1_brca2_vep:
    input:
        "data/work/{lib}/{sample}/deepvariant/brca1_brca2.output.vcf.gz"
    output:
        vcf="data/work/{lib}/{sample}/deepvariant/brca1_brca2.output.vep.vcf.gz"
    params:
        in_vcf=temp('data/work/{lib}/{sample}/deepvariant/brca1_brca2.output.vcf'),
        out_vcf='data/work/{lib}/{sample}/deepvariant/brca1_brca2.output.vep.vcf',
        assembly=config['reference']['key'],
        fa=config['reference']['fasta'],
        splice_snv=config['resources']['splice_snv'],
        splice_indel=config['resources']['splice_indel'],
        gnomAD=config['resources']['gnomAD'],
        clinvar=config['resources']['clinvar'],
        revel=config['resources']['revel'],
        loftee='$HOME/.vep/Plugins/loftee',#check
        utr=config['resources']['utrannotator']
    shell:
        """
        bcftools view -i '%FILTER=\"PASS\"' {input} | bcftools norm -m-both | bcftools norm -f {params.fa} -O v -o {params.in_vcf}

        #if {{ conda env list | grep 'vep'; }} >/dev/null 2>&1; then source activate vep; fi

        vep -i {params.in_vcf} -o {params.out_vcf} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf \
        --everything \
        --canonical \
        --assembly {params.assembly} \
        --species homo_sapiens \
        --fasta {params.fa} \
        --vcf_info_field ANN \
        --plugin NMD \
        --plugin Downstream \
        --plugin REVEL,{params.revel} \
        --plugin SpliceAI,snv={params.splice_snv},indel={params.splice_indel} \
        --plugin gnomADc,{params.gnomAD} \
        --plugin LoF,loftee_path:{params.loftee} \
        --plugin UTRannotator,{params.utr} \
        --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN

        bgzip {params.out_vcf} && tabix -fp vcf {output.vcf}

        """

rule brca1_brca2_vep_report:
    input:
        vcf="data/work/{lib}/{sample}/deepvariant/brca1_brca2.output.vep.vcf.gz"
    output:
        csv="data/work/{lib}/{sample}/deepvariant/brca1_brca2.output.vep.report.csv"
    shell:
        """
        python vep_vcf_parser.py -i {input.vcf} -o {output.csv} --mode single,{wildcards.sample} everything
        """

rule concat_reports:
    input:
        expand("data/work/{lib}/{sample}/deepvariant/brca1_brca2.output.vep.report.csv",lib=config['resources']['targets_key'],sample=SAMPLES)
    output:
        "{project}.deepvariant.brca1_brca2.output.vep.report.csv"
    shell:
        """
        awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} > {output}
        """
