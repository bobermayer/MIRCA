import os
import sys
import numpy as np
import pandas as pd
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-i','--inf',dest='inf',help="""input file""")
parser.add_option('-o','--out',dest='out',help="""output file""")
parser.add_option('-t','--type',dest='type',default='bam',help="""file type (bam|fastq|fasta) [bam]""")
parser.add_option('','--target_GC',dest='target_GC',help="""target GC distribution (output of get_read_GC_content.py)""")
parser.add_option('','--input_GC',dest='input_GC',help="""input GC distribution (output of get_read_GC_content.py)""")
parser.add_option('-n','--nbins',dest='nbins',default=10,type=int,help="number of bins for histogram")
parser.add_option('-s','--seed',dest='seed',default=0,type=int,help="""random seed [0]""")
parser.add_option('','--genome',help="""genome (indexed fasta) to use genomic (not read) sequence (only for mapped reads)""")

options,args=parser.parse_args()

np.random.seed(options.seed)

input_GC=pd.read_csv(options.input_GC,sep='\t',header=0,index_col=0,comment='#').squeeze()
target_GC=pd.read_csv(options.target_GC,sep='\t',header=0,index_col=0,comment='#').squeeze()

ratio=(target_GC/target_GC.sum())/(input_GC/input_GC.sum())
# normalize by biggest ratio (but consider only bins with appreciable amount of reads to reduce outlier influence)
ratio=ratio/ratio[input_GC > .1*input_GC.sum()/options.nbins].max()

bins=ratio.index

if options.type=='bam':
    import pysam
    inf=pysam.Samfile(options.inf)
    outf=pysam.Samfile(options.out,'wb',template=inf)
    if options.genome is not None:
        genome=pysam.FastaFile(options.genome)
else:
    if options.genome is not None:
        raise Exception("cannot get genomic read sequence for fastq!")
    import gzip
    from Bio import SeqIO
    inf=SeqIO.parse(gzip.open(options.inf,'rb') if options.inf.endswith('.gz') else options.inf,options.type)
    outf=gzip.open(options.out,'wb') if options.out.endswith('.gz') else open(options.out,'w')


GC_old=[]
GC_new=[]
for n,read in enumerate(inf):
    if n%1000==0:
        print >> sys.stderr, 'processed {0}k reads\r'.format(n/1000),
    if options.genome is not None:
        chunks=np.concatenate([[read.reference_start],read.reference_start+np.where(np.diff(read.positions) > 1)[0],[read.reference_end]])
        seq=''.join([genome.fetch(reference=read.reference_name,start=chunks[i],end=chunks[i+1]) for i in range(len(chunks)-1)]).upper()
    else:
        seq=str(read.seq).upper()
    gc=(seq.count('G')+seq.count('C'))/float(len(seq))
    i=np.searchsorted(bins,gc)
    GC_old.append(gc)
    if np.random.rand() < ratio[bins[max(i-1,0)]]:
        GC_new.append(gc)
        if options.type=='bam':
            outf.write(read)
        else:
            SeqIO.write(read,outf,options.type)

print >> sys.stderr, 'processed {0} reads'.format(n)

counts,bins=np.histogram(GC_new,bins=np.linspace(0,1,options.nbins+1))

print '# {0}: {1} of {2} reads sampled ({3:.1f}%)'.format(options.inf,len(GC_new),len(GC_old),100*len(GC_new)/float(len(GC_old)))
print '# new GC: mu={0:.3f} sigma={1:.3f}'.format(100*np.mean(GC_new),100*np.std(GC_new))
print 'GC\tfrac_reads'
print '\n'.join('{0}\t{1:.6f}'.format(b,c/float(len(GC_new))) for c,b in zip(counts,bins))
