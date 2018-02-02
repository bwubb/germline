#!/usr/bin/env python


import sys
import csv
import argparse
import re
from collections import defaultdict


class GeneError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class TranscriptError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class FileHeaderError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class CleverName(object):
    #getattr(object, name[, default])
    def __init__(self,row_dict,pTX):#Stop this rHeader,header bs. Find a way to trim the positions you need
        #print row_dict
        for attr in ['Chr','Start','End','Ref','Alt']:
            #print '{0}:{1}:{2} {3}>{4}'.format(*[row_dict[x] for x in ['Chr','Start','End','Ref','Alt']])
            setattr(self,attr,row_dict[attr])
        else:
            self.GenomicLocation='{0}:{1}-{2}'.format(self.Chr,self.Start,self.End)

        fields={'GenomicRegion':[],'IntergenicDist':[],'Gene':[],'Transcript':[],'Exon':[],'NTChange':[],'AAChange':[]}
        
        if row_dict['Func_refGene'].startswith('ncRNA'):#Dropping ncRNA for now
            pass
        
        elif row_dict['Func_refGene']=='exonic':
            eg=row_dict['Gene_refGene'].split(',')
            aa=row_dict['AAChange_refGene'].split(',')
            fields=self._exonic_genes(fields,eg,aa,pTX)
            if len(fields['Gene'])!=0:
                setattr(self,'ExonicFunc_refGene',row_dict['ExonicFunc_refGene'])
            
        elif row_dict['Func_refGene']=='exonic;splicing':
            eg=row_dict['Gene_refGene'].split(';')[0].split(',')#I think this should be whatever is left after gene filter
            aa=row_dict['AAChange_refGene'].split(',')
            fields=self._exonic_genes(fields,eg,aa,pTX)
            if len(fields['Gene'])!=0:
                setattr(self,'ExonicFunc_refGene',row_dict['ExonicFunc_refGene'])
            
            sp=row_dict['Gene_refGene'].split(';')[1].split(',')
            gd=row_dict['GeneDetail_refGene'].split(',')
            fields=self._splicing_genes(fields,sp,gd,pTX)
            
        elif row_dict['Func_refGene']=='splicing':
            sp=row_dict['Gene_refGene'].split(',')
            gd=row_dict['GeneDetail_refGene'].split(',')
            fields=self._splicing_genes(fields,sp,gd,pTX)
            
        #ncRNA_splicing and ncRNA_exonic;splicing behave almost like splicing
        #Func.refGene    Gene.refGene    GeneDetail.refGene
        #ncRNA_exonic;splicing    SENP3-EIF4A1;SENP3    NM_015670:exon7:c.1305+1A>-,NM_015670:exon8:c.1306-1A>-
        #ncRNA_exonic;splicing    ASB16-AS1;ASB16-AS1    NR_049730:exon3:c.252-1G>C
        #ncRNA_exonic;splicing    THAP7-AS1;TUBA3FP    NR_003608:exon3:c.714-2A>G
        
        #ncRNA_splicing    CLCA3P    NR_024604:exon4:c.575+2T>C

        elif row_dict['Func_refGene'] in ['UTR3','UTR5']:
            utr=row_dict['Func_refGene']
            ug=row_dict['Gene_refGene'].split(',')
            gd=row_dict['GeneDetail_refGene'].split(',')
            #Above will fail on DNAJC11,THAP3 NM_018198:c.*1161C>T;NM_138350:c.*426G>A
            #Need to switch to re.split with multiple things. Would that work and be best for all cases?
            fields=self._UTR_genes(fields,ug,gd,utr,pTX)
            
        elif row_dict['Func_refGene'] in ['UTR3;UTR5','UTR5;UTR3']:
            func=row_dict['Func_refGene'].split(';')
            for f in range(len(func)):
                ug=row_dict['Gene_refGene'].split(';')[f].split(',')
                gd=row_dict['GeneDetail_refGene'].split(';')[f].split(',')
                fields=self._UTR_genes(fields,ug,gd,func[f],pTX)
            
        elif row_dict['Func_refGene']=='intergenic':
            ig=row_dict['Gene_refGene'].split(',')
            gd=[d.lstrip('dist=') for d in row_dict['GeneDetail_refGene'].split(';')]
            fields=self._intergenic_genes(fields,ig,gd,pTX)
            
        elif row_dict['Func_refGene']=='upstream;downstream':#? downsteam;upstream?
            func=row_dict['Func_refGene'].split(';')
            for f in range(len(func)):
                og=row_dict['Gene_refGene'].split(';')[f].split(',')
                fields=self._other_genes(fields,og,func[f],pTX)
        
        else:#intronic,ncRNA_exonic,ncRNA_intronic
            func=row_dict['Func_refGene']
            og=row_dict['Gene_refGene'].split(',')
            fields=self._other_genes(fields,og,func,pTX)
        
        #if transcript is empty, fill it with pTX
        if len(fields['Gene'])==0:
            raise GeneError('GeneError: No gene returned - {0}:{1}:{2} {3}>{4} '.format(*[getattr(self,x) for x in ['Chr','Start','End','Ref','Alt']]) + row_dict['Gene_refGene'])
        else:
            for k,v in fields.items():
                setattr(self,k+'_refGene',';'.join(v))
        
    def _exonic_genes(self,fields,eg,aa,pTX):
        g,a,n,e,t=([] for x in range(5))
        for i in eg:
            V=None
            for j in aa:
                if pTX[i] in j:
                    V=j.split(':')
            if V:
                assert len(V)>3,'AssertionError: _exonic_genes - gene_info = {0}'.format(V)#Nonframeshift substitutions do not have amino annotation
                try:
                    a.append(V[4])
                except IndexError:
                    a.append('')
                n.append(V[3])
                e.append(V[2][4:])
                t.append(V[1])
                g.append(V[0])
        if len(g)>0:
            for k,v in zip(['Gene','Transcript','Exon','NTChange','AAChange'],[g,t,e,n,a]):
                fields[k]+=[','.join(v)]
            fields['GenomicRegion']+=['exonic']
        return fields
        
    def _splicing_genes(self,fields,sp,gd,pTX):
        g,n,e,t=([] for x in range(4))
        for i in sp:
            V=None
            for j in gd:
                if pTX[i] in j:
                    V=j.split(':')
            if V:
                assert len(V)==3,'AssertionError: _splicing_genes - gene_info = {0}'.format(V)
                n.append(V[2])
                e.append(V[1])
                t.append(V[0])
                g.append(i)
        if len(g)>0:
            for k,v in zip(['Gene','Transcript','Exon','NTChange'],[g,t,e,n]):
                fields[k]+=[','.join(v)]
            fields['GenomicRegion']+=['splicing']
        return fields
    
    def _UTR_genes(self,fields,ug,gd,utr,pTX):
        g,n,t=([] for x in range(3))
        for i in ug:
            V=None
            for j in gd:
                if pTX[i] in j:
                    V=j.split(':')
            if V:
                assert len(V)==2,'AssertionError: _UTR_genes - gene_info = {0}'.format(V)
                n.append(V[1])
                t.append(V[0])
                g.append(i)
        if len(g)>0:
            for k,v in zip(['Gene','Transcript','NTChange'],[g,t,n]):
                fields[k]+=[','.join(v)]
            fields['GenomicRegion']+=[utr]
        return fields
    
    def _intergenic_genes(self,fields,ig,gd,pTX):
        d,g,t=([] for x in range(3))
        for n,i in enumerate(ig):
            if pTX[i]:
                g.append(i)
                d.append(':'.join([i,gd[n]]))
                t.append(pTX[i])
        if len(g)>0:
            for k,v in zip(['IntergenicDist','Gene','Transcript'],[d,g,t]):
                fields[k]+=[','.join(v)]
            fields['GenomicRegion']+=['intergenic']
        return fields
    
    def _other_genes(self,fields,og,func,pTX):
        g,t=([] for x in range(2))
        for i in og:
            if pTX.get(i,False):#####Look at this notation! Otherwise you can get a key error from shit like ASB16-AS1;ASB16-AS1
                g.append(i)
                t.append(pTX[i])
        if len(g)>0:
            for k,v in zip(['Gene','Transcript'],[g,t]):
                fields[k]+=[','.join(v)]
            fields['GenomicRegion']+=[func]
        return fields
    
    def _region_info(self,i):
        if i['genomicSuperDups']:
            #genomicSuperDups.Score
            #genomicSuperDups.Name
            v=i['genomicSuperDups'].split(';')
            if len(v)==2:
                setattr(self,'genomicSuperDups_Score',v[0].split('=')[1])
                setattr(self,'genomicSuperDups_Name',v[1].split('=')[1])
    
    def _filter_info(self,i,header):
        for h in i.keys():
            if h in header and not hasattr(self,h):
                setattr(self,h,i[h])
            elif h.startswith('snp'):#snp138,snp138NonFlagged
                ver=re.match('snp(\d{3})',h).group(1)
                if h.endswith('NonFlagged'):
                    setattr(self,'dbSNP_NonFlagged',i[h])
                else:
                    setattr(self,'dbSNP',i[h])
                    setattr(self,'dbSNP_ver',ver)#No overwrite
            elif h.startswith('cosmic'):
                ver=re.match('cosmic(\d+)',h).group(1)
                if h.endswith('wgs'):
                    setattr(self,'COSMICwgs',i[h])#switch to . or _ notation
                else:
                    setattr(self,'COSMIC',i[h])
                    setattr(self,'COSMIC_ver',ver)#No overwrite
            elif h.startswith('CL'):
                k='ClinVar_'+h
                setattr(self,k,i[h])
            elif h=='nci60':
                setattr(self,'NCI-60',i[h])
            elif h=='Interpro_domain':
                setattr(self,'dbNSFP_InterPro_domain',i[h])
            elif h=='dbscSNV_ADA_SCORE':
                setattr(self,'dbscSNV_ada_score',i[h])
            elif h=='dbscSNV_RF_SCORE':
                setattr(self,'dbscSNV_rf_score',i[h])
                #Now this only works with using getattr
            else:
                pass
    
    def _coverage_info(self,i):
        DP=0
        if i.get('GT'):
            i['GT']=i['GT'].replace('/','|')
            try:
                self.Zyg={'.|.':'','0':'hom_ref','1':'hom_alt','.':'unknown','1|1':'hom_alt','0|1':'het','1|0':'het','0|0':'hom_ref'}[i['GT']]
            except KeyError:
                self.Zyg=i['GT']+'_alt'
        if i.get('DP'):
            DP=i['DP']
            setattr(self,'Depth',i['DP'])
        if i.get('AD'):
            AD=i['AD'].split(',')[1:]
            AD=[a for a in AD if int(a)!=0]
            setattr(self,'ALT_AlleleDepth',','.join(AD))
        elif i.get('AO'):
            AD=i['AO'].split(',')
            AD=[a for a in AD if int(a)!=0]
            setattr(self,'ALT_AlleleDepth',','.join(AD))
        else:
            print('No Alt allele count found',i)
        if getattr(self,'Depth','')=='':
            try:
                DP=str(sum([int(a) for a in i['AD'].split(',')]))
                setattr(self,'Depth',DP)
            except (ZeroDivisionError,TypeError,ValueError) as e:
                pass
        try:
            setattr(self,'ALT_AlleleFrac',','.join(['{0:.3f}'.format(float(ad)/float(DP)) for ad in AD]))
        except (ZeroDivisionError,TypeError,ValueError) as e:
            print('{0}: DP ={1} AD ={2}  {3}'.format(e,DP,AD,i))
    
    def _targeted_status(self,targets=None):
        if any(self._in_range(s,p,e) for s,e in targets[self.Chr] for p in [int(self.Start),int(self.End)]):
            self.Targeted_status='{0};OnTarget'.format(self.SeqLibrary)
        else:
            self.Targeted_status='{0};OffTarget'.format(self.SeqLibrary)
    
    def _sample_info(self,info):#"meta" data
        for k,v in info.items():
            setattr(self,k,v)
    
    def _matched_variant(self,i):
        if i.get('DP'):
            DP=i['DP']
            setattr(self,'MATCHED_Depth',i['DP'])
        if i.get('AD'):
            AD=i['AD'].split(',')[1:]
            AD=[a for a in AD if int(a)!=0]
            setattr(self,'MATCHED_ALT_AlleleDepth',','.join(AD))
        elif i.get('AO'):
            AD=i['AO'].split(',')
            AD=[a for a in AD if int(a)!=0]
            setattr(self,'MATCHED_ALT_AlleleDepth',','.join(AD))
        else:
            print('No Alt allele count found',i)
        setattr(self,'MATCHED_ALT_AlleleFrac',','.join(['{0:.3f}'.format(float(ad)/float(DP)) for ad in AD]))
        try:
            setattr(self,'MATCHED_ALT_AlleleFrac',','.join(['{0:.3f}'.format(float(ad)/float(DP)) for ad in AD]))
        except (ZeroDivisionError,TypeError,ValueError) as e:
            print('{0}: DP ={1} AD ={2}  {3}'.format(e,DP,AD,i))
    
    def __repr__(self):
        x=defaultdict(str)
        for i in ['Chr','Start','End','Gene_refGene','Transcript_refGene']:
            x[i]=getattr(self,i,'')
        return str(x)
    
    def _in_range(self,s,p,e):
        return s-100<=p<=e+100
    
    def out(self,header):
        return [getattr(self,h,str()) for h in header]

#####

def get_principleTX():
    pTX={}
    with open('/home/bwubb/resources/preferred_transcripts.20170505.txt','r') as file:
        for line in file:
            try:
                k,v=line.rstrip().split('\t')
                pTX[k]=v
            except IndexError:
                pass
    return pTX

def get_header():
    header='SampleID ProjectID SeqLibrary VariantCaller Chr Start End Ref Alt GenomicLocation'.split(' ')
    header+='GenomicRegion_refGene ExonicFunc_refGene Gene_refGene Transcript_refGene Exon_refGene NTChange_refGene AAChange_refGene'.split(' ')
    header+='dbNSFP_InterPro_domain cytoband gwasCatalog genomicSuperDups_Name genomicSuperDups_Score dbscSNV_ada_score dbscSNV_rf_score'.split(' ')
    header+='SIFT_score SIFT_pred Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred MutationAssessor_score MutationAssessor_pred'.split(' ')
    header+='FATHMM_score FATHMM_pred PROVEAN_score PROVEAN_pred VEST3_score CADD_raw CADD_phred DANN_score fathmm-MKL_coding_score fathmm-MKL_coding_pred MetaSVM_score MetaSVM_pred MetaLR_score MetaLR_pred integrated_fitCons_score integrated_confidence_value GERP++_RS phyloP7way_vertebrate phyloP20way_mammalian phastCons7way_vertebrate phastCons20way_mammalian SiPhy_29way_logOdds'.split(' ')
    header+='dbSNP dbSNP_NonFlagged dbSNP_ver COSMIC COSMIC_ver PopFreqMax'.split(' ')
    header+='1000G_ALL 1000G_AFR 1000G_AMR 1000G_EAS 1000G_EUR 1000G_SAS'.split(' ')
    header+='ExAC_ALL ExAC_AFR ExAC_AMR ExAC_EAS ExAC_FIN ExAC_NFE ExAC_OTH ExAC_SAS'.split(' ')
    header+='ExAC_nontcga_ALL ExAC_nontcga_AFR ExAC_nontcga_AMR ExAC_nontcga_EAS ExAC_nontcga_FIN ExAC_nontcga_NFE ExAC_nontcga_OTH ExAC_nontcga_SAS'.split(' ')
    header+='ESP6500siv2_ALL ESP6500siv2_AA ESP6500siv2_EA CG46 NCI-60'.split(' ')
    header+='ClinVar_CLINSIG ClinVar_CLNDBN ClinVar_CLNACC'.split(' ')
    header+='Zyg Depth ALT_AlleleDepth ALT_AlleleFrac'.split(' ')
    return header

def summarize_runner(reader,writer,header,fileinfo,pTX):
    writer.writerow(header)
    raw_header=next(reader)#header line
    raw_header[-1]='QUAL'
    raw_header+=['FILTER','INFO','FORMAT','SAMPLE']#header has less columns than rows; carry over from vcf
    raw_header=[i.replace('.','_') for i in raw_header]
    for row in reader:
        row=fix_blanks(row)
        row_dict=dict(zip(raw_header,row))
        try:
            genes=list(set(row_dict['Gene_refGene'].replace(';',',').split(',')))#'GeneA;GeneB,GeneA' > GeneA,GeneB,GeneA > [GeneA,GeneB]
            #can you figure out a clever way to makes tx,tx;tx?
            if not all(gene in pTX.keys() for gene in genes): #This part is good, need exception for rs ids that are in the design. They could be noncoding rna
                raise TranscriptError("{0} {1} {2} {3} {4} {5} {6}".format(row_dict['Chr'],row_dict['Start'],row_dict['End'],row_dict['Ref'],row_dict['Alt'],row_dict['Gene_refGene'],row_dict['AAChange_refGene']))
            Data=CleverName(row_dict,pTX)
            Data._filter_info(row_dict,header)
            cov_dict=dict(zip(row_dict['FORMAT'].split(':'),row_dict['SAMPLE'].split(':')))
            Data._coverage_info(cov_dict)
            Data._sample_info(fileinfo)
        except (GeneError,TranscriptError,AssertionError,UnboundLocalError) as e:
            print(e)
            continue
        writer.writerow(Data.out(header))#need meta header too.

def pull(x,field,value):
    i=defaultdict(bool)
    for h in field:
        if any(h.startswith(y) for y in x):
            i[h]=value[field.index(h)]
    return i

def fix_blanks(x):
    for i,j in enumerate(x):
        if j=='.':
            x[i]=''
    return x

def get_args():
    p=argparse.ArgumentParser()
    p.add_argument('--input_fp',help='File path for input')
    p.add_argument('--output_fp',help='File path for output')
    p.add_argument('--SM',help='SampleID')
    p.add_argument('--PR',help='ProjectID')
    p.add_argument('--LB',help='SeqLibrary')
    p.add_argument('--VC',help='VariantCaller')
    return p.parse_args()

def main(argv=None):
    csv.field_size_limit(sys.maxsize)
    pTX=get_principleTX()
    header=get_header()
    args=get_args()
    fileinfo={'SampleID':args.SM,'ProjectID':args.PR,'SeqLibrary':args.LB,'VariantCaller':args.VC}
    with open(args.input_fp,'r') as input, open(args.output_fp,'w') as output:
            reader=csv.reader(input,delimiter='\t')#Need to convert to DictReader
            writer=csv.writer(output,delimiter='\t')#Need to convert to DictWriter
            summarize_runner(reader,writer,header,fileinfo,pTX)

if __name__ == '__main__':
    main()
