#!/usr/bin/env python

import os
import argparse
import errno

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno==errno.EEXIST and os.path.isdir(path):
			pass

def protocols_operations():
	params={}
	params['protocol']=[]
	params['operation']=[]
	with open('','rb') as file:
		for line in file:
			if not line.startswith('#'):
				p,o=line.rstrip().split('\t')
				params['protocol'].append(p)
				params['operation'].append(o)
				params[p]=o
	return params

def write_script(sample,vcf,intervals):
	out_vcf=os.path.abspath('data/final/{0}/{0}.{1}.gatk-haplotype.vcf'.format(sample,os.path.basename(os.getcwd())))
	out_a=os.path.abspath('data/annotation/{0}/{0}.{1}.gatk-haplotype'.format(sample,os.path.basename(os.getcwd())))
	#out='{0}.haplotype'.format(sample)
	with open('data/config/{0}/gatk-annotate.sh'.format(sample),'wb') as script:
		script.write('#!/bin/bash\n\n')
		script.write('DATE=`date "+%Y%m%d"`\n')
		script.write('ref="/home/bwubb/b37_genomes/human_g1k_v37.fasta"\n')
		script.write('GATK="/home/bwubb/GenomeAnalysisTK-3.4-46"\n')
		script.write('humandb="/home/bwubb/humandb"\n')
		#script.write('zcat {0}.vcf.gz > {0}.vcf\n\n'.format(out))
		script.write('java -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \\\n')
		script.write(' -T SelectVariants \\\n')
		script.write(' -R $ref \\\n')
		script.write(' -V {0} \\\n'.format(vcf))
		if intervals:
			script.write(' -L {0} \\\n'.format(intervals))
		script.write(' -sn {0} \\\n'.format(sample))
		script.write(' -env -ef \\\n')
		script.write(' -o {0}\n\n'.format(out_vcf))
		
		script.write('convert2annovar.pl -format vcf4 -includeinfo {0} | awk -F"\\t" '.format(out_vcf)+"'BEGIN { OFS=FS } { print $1,$2,$3,$4,$5,$11,$12,$13,$14,$15 }' > "+'{0}.avinput\n'.format(out_a))
		script.write('table_annovar.pl {0}.avinput $humandb -build hg19 --outfile {0} -protocol refGene,cytoband,gwasCatalog,genomicSuperDups,dbnsfp31a_interpro,dbscsnv11,dbnsfp30a,snp138,snp138NonFlagged,cosmic70,popfreq_all_20150413,exac03nontcga,nci60,clinvar_20160302 -operation g,r,r,r,f,f,f,f,f,f,f,f,f,f -otherinfo -remove\n\n'.format(out_a))

#Ideally I would write this script to be executed directly instead of bash

def read_input(infile):
	with open(infile,'rb') as input:
		samples=input.read().splitlines()
	for sample in samples:
		vcf=None
		for path,dirs,files in os.walk('data/final/{0}'.format(sample)):
			for file in files:
				if file=='{0}-gatk-haplotype.vcf.gz'.format(sample):
					vcf=os.path.abspath(os.path.join(path,file))
					print vcf
					write_script(sample,vcf)
		if not vcf:
			print 'ERROR: Could not locate {0}-gatk-haplotype.vcf.gz'.format(sample)

def get_args():
	p=argparse.ArgumentParser()
	p.add_argument('-i','--infile',help='List of samples')
	p.add_argument('-V','--vcf',help='VCF file')
	p.add_argument('-L','--intervals',help='')
	
	#p.add_argument('--build',help='buildver')
	#p.add_argument('--gatk-version')
	return p.parse_args()

def main(argv=None):
	args=get_args()
	for a,b in vars(args).items():
		print '{}: {}'.format(a,b)
	
	with open(args.infile,'rb') as file:
		samples=file.read().splitlines()
	mkdir_p('data/annotation')
	for sample in samples:
		mkdir_p('data/annotation/{0}'.format(sample))
		write_script(sample,os.path.abspath(args.vcf),args.intervals)
	#read_input(args.infile)

if __name__ == '__main__':
	main()