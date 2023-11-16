import os

with open(config.get('project',{}).get('sample_list','sample.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config.get('project',{}).get('bam_table','bam.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

def map_preprocess(wildcards):
    return {'bam':BAMS[wildcards.sample],'bcf':f'data/work/{wildcards.lib}/{wildcards.sample}/manta-wes/{wildcards.candidates}.bcf','aln':f'data/work/{wildcards.lib}/{wildcards.sample}/varlociraptor/alignment-properties.json'}

def candidate_header(wildcards):
    #H={}
    return f"$HOME/resources/Vcf_files/simplified-header-w_contig.{config['reference']['key'].lower()}.vcf"

wildcard_constraints:
    candidates="\w+"
#    caller="^manta-\w+",
#    candidates="^candidate\w+"

rule manta_fdr:
    input: expand("data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/{candidates}.local-fdr.bcf",lib=config['resources']['targets_key'],sample=SAMPLES,candidates=["candidateSmallIndels","candidateSV"])

rule manta_annotation:
    input: expand("data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/candidateSmallIndels.local-fdr.vep.vcf.gz",lib=config['resources']['targets_key'],sample=SAMPLES),
        expand("data/work/{lib}/{sample}/annotsv/manta-wes/candidateSV.local-fdr.annotated.tsv",lib=config['resources']['targets_key'],sample=SAMPLES,candidates=["candidateSmallIndels","candidateSV"])

#Manta has wgs or wes flavours. For now I still think it is important to distinguish between the two.
rule manta_candidate_tsv:
    input:
        #"data/work/{lib}/{sample}/{caller}/results/variants/candidateSmallIndels.vcf.gz",
        #"data/work/{lib}/{sample}/{caller}/results/variants/candidateSV.vcf.gz"
        "data/work/{lib}/{sample}/manta-wes/results/variants/{candidates}.vcf.gz"
    output:
        "data/work/{lib}/{sample}/manta-wes/{candidates}.tsv"
    shell:
        """
        bcftools concat -a {input} | bcftools query -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t%FILTER\t.\n' > {output}
        """

rule manta_candidate_bcf:
    input:
        "data/work/{lib}/{sample}/manta-wes/{candidates}.tsv"
    output:
        "data/work/{lib}/{sample}/manta-wes/{candidates}.bcf"
    params:
        tsv=temp("data/work/{lib}/{sample}/manta-wes/tmp.{candidates}.tsv"),
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
        "data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/{candidates}.observations.bcf"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        varlociraptor preprocess variants {params.ref} --alignment-properties {input.aln} --bam {input.bam} --candidates {input.bcf} > {output}
        """

rule varlociraptor_scenario:
    input:
        yml='germline_scenario.yml',
        bcf="data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/{candidates}.observations.bcf"
    output:
        "data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/{candidates}.bcf"
    shell:
        """
        varlociraptor call variants generic --scenario {input.yml} --obs germline={input.bcf} > {output}
        """

#this leaves behind bcf index files.
rule varlociraptor_local_fdr:
    input:
        "data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/{candidates}.bcf"
    output:
        "data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/{candidates}.local-fdr.bcf"
    params:
        fdr="0.05",
        het=temp("data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/{candidates}.local-fdr.het.bcf"),
        hom=temp("data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/{candidates}.local-fdr.hom.bcf")
    shell:
        """
        varlociraptor filter-calls control-fdr --local {input} --events HET --fdr {params.fdr} > {params.het}
        bcftools index {params.het}
        varlociraptor filter-calls control-fdr --local {input} --events HOM --fdr {params.fdr} > {params.hom}
        bcftools index {params.hom}
        bcftools concat -a {params.het} {params.hom} | bcftools sort -O b -o {output}
        bcftools index {output}
        """

rule annotsv_candidateSV:
    input:
        "data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/candidateSV.local-fdr.bcf"
    output:
        "data/work/{lib}/{sample}/annotsv/manta-wes/candidateSV.local-fdr.annotated.tsv",
        "data/work/{lib}/{sample}/annotsv/manta-wes/candidateSV.local-fdr.unannotated.tsv"
    params:
        in_vcf=temp("data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/candidateSV.local-fdr.vcf"),
        out_dir="data/work/{lib}/{sample}/annotsv/manta-wes/",
        build=config['reference']['key'],
        tx="ENSEMBL" #config['analysis']['tx']
    shell:
        """
        bcftools view -O v -o {params.in_vcf} {input}
        AnnotSV -SVinputFile {params.in_vcf} -outputDir {params.out_dir} -annotationMode split -genomeBuild {params.build}
        """


rule vep_candidateSmallIndels:
    input:
        "data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/candidateSmallIndels.local-fdr.bcf"
    output:
        "data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/candidateSmallIndels.local-fdr.vep.vcf.gz"
    params:
        in_vcf=temp('data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/candidateSmallIndels.local-fdr.vcf'),
        out_vcf=temp('data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/candidateSmallIndels.local-fdr.vep.vcf'),
        assembly=config['reference']['key'],
        fa=config['reference']['fasta'],
        splice_snv=config['resources']['splice_snv'],
        splice_indel=config['resources']['splice_indel'],
        gnomAD=config['resources']['gnomAD'],
        clinvar=config['resources']['clinvar'],
        revel=config['resources']['revel'],
        loftee='$HOME/.vep/Plugins/loftee',#check
        utr=config['resources']['utr'],
        ref_fa="data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/reference.fa",
        mut_fa="data/work/{lib}/{sample}/varlociraptor/germline_scenario/manta-wes/mutated.fa"
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

        bgzip -c {params.out_vcf} > {output}
        tabix -fp vcf {output}
        """
