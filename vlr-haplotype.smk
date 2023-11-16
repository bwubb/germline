import os
#include ./lancet.smk
#include ./mutect2.smk
#include ./strelka2.smk
#include ./vardictjava.smk
#include ./varscan2.smk

##INIT
with open(config.get('project',{}).get('sample_list','sample.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config.get('project',{}).get('bam_table','bam.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

CHR=[f'chr{x}' for x in range(1,23)]+['chrX']

##PYTHON
def map_vcf(wildcards):
    V={}
    return V[wildcards.caller]

def map_preprocess(wildcards):
    return {'bam':BAMS[wildcards.sample],'bcf':f'data/work/{wildcards.lib}/{wildcards.sample}/varlociraptor/candidates.{wildcards.caller}.bcf','aln':f'data/work/{wildcards.lib}/{wildcards.sample}/varlociraptor/alignment-properties.json'}

def candidate_header(wildcards):
    #H={}
    return "$HOME/resources/Vcf_files/simplified-header-w_contig.grch38.vcf"

def genome_size(wildcards):
    G={'GRCh38':'3.3e9','S04380110':'5.0e7','S07604715':'6.6e7','S31285117':'4.9e7','xgen-exome-research-panel-targets-grch37':'3.9e7'}
    return G[config['resources']['targets_key']]

##TARGET RULES
wildcard_constraints:
    chr="chr[0-9XY]+",
    caller="[a-z]+",
    scenario="germline_scenario"

rule collect_haplotype:
    input:
        expand("data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.bcf",lib=config['resources']['targets_key'],sample=SAMPLES,scenario="germline_scenario",caller="gatk")

rule collect_haplotype_vep:
    input:
        expand("data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.vep.report.csv",lib=config['resources']['targets_key'],sample=SAMPLES,scenario="germline_scenario",caller="gatk")

rule collect_manta:
    input:
        expand("data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.bcf",lib=config['resources']['targets_key'],sample=SAMPLES,scenario="germline_scenario",caller="manta")

rule collect_manta_vep:
    input:
        expand("data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.vep.report.csv",lib=config['resources']['targets_key'],sample=SAMPLES,scenario="germline_scenario",caller="manta")


##SNAKEMAKE
rule gatk_candidate_tsv_byChr:
    input:
        expand("data/work/{{lib}}/{project}/gatk/haplotype.{{chr}}.vcf.gz",project=config['project']['name'])
    output:
        "data/work/{lib}/{sample}/gatk/candidates.{chr}.tsv"
    shell:
        """
        bcftools view -a -s {wildcards.sample} -e 'ALT~\"*\"' {input} | bcftools query -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t%FILTER\t.\n' > {output}
        """
    #removed -i 'FILTER=\"PASS\"' from query

rule gatk_candidate_bcf_byChr:
    input:
        expand("data/work/{{lib}}/{{sample}}/gatk/candidates.{chr}.tsv",chr=CHR)
    output:
        "data/work/{lib}/{sample}/varlociraptor/candidates.gatk.bcf"
    params:
        tsv=temp("data/work/{lib}/{sample}/varlociraptor/tmp.candidates.gatk.tsv"),
        header=candidate_header
    shell:
        """
        cat {input} | sort -u > {params.tsv}
        cat {params.header} {params.tsv} | bcftools sort -O b -o {output}
        bcftools index {output}
        """

rule manta_candidate_tsv:
    input:
        "data/work/{lib}/{sample}/manta-wgs/results/variants/candidateSmallIndels.vcf.gz",
        "data/work/{lib}/{sample}/manta-wgs/results/variants/candidateSV.vcf.gz"
    output:
        "data/work/{lib}/{sample}/manta-wgs/candidates.tsv"
    shell:
        """
        bcftools concat -a {input} | bcftools query -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t%FILTER\t.\n' > {output}
        """

rule manta_candidate_bcf:
    input:
        "data/work/{lib}/{sample}/manta-wgs/candidates.tsv"
    output:
        "data/work/{lib}/{sample}/varlociraptor/candidates.manta.bcf"
    params:
        tsv=temp("data/work/{lib}/{sample}/varlociraptor/tmp.candidates.manta.tsv"),
        header=candidate_header
    shell:
        """
        cat {input} | sort -u > {params.tsv}
        cat {params.header} {params.tsv} | bcftools sort -O b -o {output}
        bcftools index {output}
        """

rule varlociraptor_estimate_properties:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        "data/work/{lib}/{sample}/varlociraptor/alignment-properties.json"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        varlociraptor estimate alignment-properties {params.ref} --bam {input.bam} > {output}
        """

rule varlociraptor_preprocess_sample:
    input:
        unpack(map_preprocess)
    output:
        "data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.observations.bcf"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        varlociraptor preprocess variants {params.ref} --alignment-properties {input.aln} --bam {input.bam} --candidates {input.bcf} > {output}
        """

rule varlociraptor_scenario:
    input:
        yml='germline_scenario.yml',
        bcf="data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.observations.bcf"
    output:
        "data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.bcf"
    shell:
        """
        varlociraptor call variants generic --scenario {input.yml} --obs germline={input.bcf} > {output}
        """

#this leaves behind bcf index files.
rule varlociraptor_local_fdr:
    input:
        "data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.bcf"
    output:
        "data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.bcf"
    params:
        fdr="0.05",
        het=temp("data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.het.bcf"),
        hom=temp("data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.hom.bcf")
    shell:
        """
        varlociraptor filter-calls control-fdr --local {input} --events HET --fdr {params.fdr} > {params.het}
        bcftools index {params.het}
        varlociraptor filter-calls control-fdr --local {input} --events HOM --fdr {params.fdr} > {params.hom}
        bcftools index {params.hom}
        bcftools concat -a {params.het} {params.hom} | bcftools sort -O b -o {output}
        bcftools index {output}
        """

rule varlociraptor_vep:
    input:
        "data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.bcf"
    output:
        vcf="data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.vep.vcf.gz"
    params:
        in_vcf=temp('data/work/{lib}/{sample}/varlociraptor/{scenario}/gatk/haplotype.local-fdr.vcf'),
        out_vcf=temp('data/work/{lib}/{sample}/varlociraptor/{scenario}/gatk/haplotype.local-fdr.vep.vcf'),
        assembly=config['reference']['key'],
        fa=config['reference']['fasta'],
        splice_snv=config['resources']['splice_snv'],
        splice_indel=config['resources']['splice_indel'],
        gnomAD=config['resources']['gnomAD'],
        clinvar=config['resources']['clinvar'],
        revel=config['resources']['revel'],
        loftee='$HOME/.vep/Plugins/loftee',#check
        utr=config['resources']['utr'],
        ref_fa="data/work/{lib}/{sample}/varlociraptor/{scenario}/gatk/reference.fa",
        mut_fa="data/work/{lib}/{sample}/varlociraptor/{scenario}/gatk/mutated.fa"
    shell:
        """
        bcftools view -O v -o {params.in_vcf} {input}

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
        --plugin NMD \
        --plugin ProteinSeqs,{params.ref_fa},{params.mut_fa} \
        --plugin Downstream \
        --plugin REVEL,{params.revel} \
        --plugin SpliceAI,snv={params.splice_snv},indel={params.splice_indel} \
        --plugin gnomADc,{params.gnomAD} \
        --plugin UTRannotator,{params.utr} \
        --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN

        bgzip -c {params.out_vcf} > {output.vcf}
        tabix -fp vcf {output.vcf}
        """
#        --plugin LoF,loftee_path:{params.loftee} \

rule varlociraptor_vep_report:
    input:
        vcf="data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.vep.vcf.gz"
    output:
        csv="data/work/{lib}/{sample}/varlociraptor/{scenario}.{caller}.local-fdr.vep.report.csv"
    shell:
        """
        python vep_vcf_parser.py -i {input.vcf} -o {output.csv} -t Germline --mode single,{wildcards.sample} everything
        """
