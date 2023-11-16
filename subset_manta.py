#AnnotSV subset

import csv
import argparse
from collections import defaultdict

#AnnotSV_ID      SV_chrom        SV_start        SV_end


p=argparse.ArgumentParser()
p.add_argument('-i','--infile',help='AnnotSV output file')
p.add_argument('-o','--outfile',help='Subset output')
p.add_argument('-b','--bed',help='Bed regions')
args=p.parse_args()

bed_regions=defaultdict(list)
with open(args.bed,'r') as file:
    reader=csv.reader(file,delimiter='\t')
    for row in reader:
        bed_regions[row[0]].append([int(row[1]),int(row[2])])

with open(args.infile,'r') as infile,open(args.outfile,'w') as outfile:
    reader=csv.DictReader(infile,delimiter='\t')
    header=reader.fieldnames
    writer=csv.DictWriter(outfile,delimiter='\t',fieldnames=header)
    writer.writeheader()
    for row in reader:
        intervals=bed_regions[row['SV_chrom']]
        if any([start<=int(row['SV_start'])<=end for start,end in intervals]):
            writer.writerow(row)
