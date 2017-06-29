#!/usr/bin/env python

import os
import errno
import argparse
import datetime

class Interval(object):
	def __init__(self,file,key):
		self.file=file
		self.key=key

class Sample(object):
	def __init__(self,name,bam):
		self.name=name
		self.bam=bam

def get_ref(build):
	v={'GRCh37':'/home/bwubb/b37_genomes/human_g1k_v37.fasta'}
	return v[build]

def get_dbsnp(build):
	v={'GRCh37':'/home/bwubb/b37_genomes/dbsnp_137.b37.vcf'}
	return v[build]

def get_GATK(gatkv):
	v={'3.4':'/home/bwubb/GenomeAnalysisTK-3.4-46','3.7':'/home/bwubb/software/GenomeAnalysisTK-3.7'}
	return v[gatkv]

def write_script(sample,bam_out,build,gatkv,lib,mem):
	out_gvcf=os.path.abspath('data/work/{0}/{0}.{1}.g.vcf'.format(sample.name,lib.key))
	with open(os.path.abspath('data/config/{0}/haplotype.sh'.format(sample.name)),'wb') as script:
		script.write('#!/bin/bash\n\n')
		script.write('DATE={0}\n'.format(datetime.date.today().strftime("%Y%m%d")))
		script.write('ref="{0}"\n'.format(get_ref(build)))
		script.write('dbsnp="{0}"\n'.format(get_dbsnp(build)))
		script.write('GATK="{0}"\n\n'.format(get_GATK(gatkv)))
		#
		script.write('java -Xmx{0}g -jar $GATK/GenomeAnalysisTK.jar \\\n'.format(mem))
		script.write(' -R $ref \\\n')
		script.write(' -T HaplotypeCaller \\\n')
		script.write(' -I {0} \\\n'.format(sample.bam))
		script.write(' --emitRefConfidence GVCF \\\n')
		script.write(' --dbsnp $dbsnp \\\n')
		if bam_out:
			script.write(' --bamOutput {0} \\\n'.format(os.path.abspath('data/final/{0}/{0}.gatk.bam'.format(sample.name))))
		if lib.file!=None:
			script.write(' -L {0} \\\n'.format(lib.file))
		script.write(' -o {0}\n\n'.format(out_gvcf))


def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno==errno.EEXIST and os.path.isdir(path):
			pass

def get_args():
	'''Parse sys.argv'''
	parser=argparse.ArgumentParser()
	parser.add_argument('-i','--infile', help='')#def is_file or not?
	parser.add_argument('--lib',help='Name associated with capture library.')
	parser.add_argument('--bam-out',action='store_true',default=False,help='Write bam of realignment')
	parser.add_argument('--build',default='GRCh37',help='What Reference version would you like?')
	parser.add_argument('--gatkv',default='3.4',help='If you want to change the version of GATK used.')
	parser.add_argument('-M','--mem',default=10,type=int,help='How much memory for the java stack. Measured in g')
	return parser.parse_args()


def main(argv=None):
	args=get_args()
	if args.lib!=None:
		#load intervals file
		with open('germline/intervals_available.txt','rb') as file:
			intv=dict(line.rstrip().split(',') for line in file)
		if args.lib not in intv:
			raise
		lib=Interval(intv[args.lib],args.lib)
	else:
		lib=Interval(None,'all-sites')
	
	with open(args.infile,'rb') as file:
		samples=file.read().splitlines()
	
	mkdir_p('data')
	mkdir_p('data/config')
	mkdir_p('data/final')
	mkdir_p('data/work')
	
	for name in samples:
		bam=False
	
		if os.path.isfile('bam_input/final/{0}/{0}.ready.bam'.format(name)):
			bam=os.path.abspath('bam_input/final/{0}/{0}.ready.bam'.format(name))
		elif os.path.isfile('bam_input/final/{0}/{0}-ready.bam'.format(name)):
			bam=os.path.abspath('bam_input/final/{0}/{0}-ready.bam'.format(name))
		else:
			print 'Could not locate', name
		
		if bam:
			mkdir_p('data/config/{0}'.format(name))
			mkdir_p('data/final/{0}'.format(name))
			mkdir_p('data/work/{0}'.format(name))
			sample=Sample(name,bam)
			write_script(sample,args.bam_out,args.build,args.gatkv,lib,args.mem)

if __name__=='__main__':
	main()
