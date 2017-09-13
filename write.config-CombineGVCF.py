

import numpy
import os
import argparse

def _subdivisions(l,n):
	assert isinstance(n,int)
	return numpy.array_split(numpy.array(l),n)

def write_input_list(input,n):
	with open('split_input.{0}.list'.format(n+1),'wb') as file:
		for i in input:
			file.write(i+'\n')

def find_gvcf(sample):
	if os.path.isfile('data/final/{0}/{0}.TGCT.gvcf.gz'.format(sample)):
		return os.path.abspath('data/final/{0}/{0}.TGCT.gvcf.gz'.format(sample))
	else:
		return False

def get_args():
	p=argparse.ArgumentParser()
	p.add_argument('-i','--infile',help='List of samples')
	p.add_argument('--num',type=int,default=1,help='number of subdivisions')
	return p.parse_args()

def main(argv=None):
	args=get_args()
	
	with open(args.infile,'rb') as file:
		samples=file.read().splitlines()
	
	gvcf_files={}
	for sample in samples:
		f=find_gvcf(sample)
		if f:
			gvcf_files[sample]=f
		else:
			print 'No GVCF for sample:'+sample
	for n,input in enumerate(_subdivisions(gvcf_files.keys(),args.num)):
		write_input_list([gvcf_files[i] for i in input],n)


if __name__=='__main__':
	main()