import os
import sys
import gzip
import scipy.stats
import numpy as np
import pysam
from Bio import SeqIO
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-i','--in',dest='inf',help="input file (fastq or bam)")
parser.add_option('-n','--nbins',dest='nbins',default=10,type=int,help="number of bins for histogram")
parser.add_option('-t','--type',dest='type',default='bam',help="""input file type (bam / fastq / fasta)""")
parser.add_option('','--nmax',dest='nmax',type=int,help="""max number of reads to consider""")
parser.add_option('','--regress_multiplicity',action='store_true',default=False,help="""regress GC content on read multiplicity (input is output of fastx_collapser)""")
parser.add_option('','--genome',help="""genome (indexed fasta) to get GC stats of genomic (not read) sequence (only for mapped reads)""")

options,args=parser.parse_args()

if options.type=='bam':
    inf=pysam.Samfile(options.inf)
    if options.genome is not None:
        genome=pysam.FastaFile(options.genome)
else:
    if options.genome is not None:
        raise Exception("cannot get genomic read sequence for fastq!")
    if options.inf.endswith('.gz'):
        inf=SeqIO.parse(gzip.open(options.inf,'rb'),options.type)
    else:
        inf=SeqIO.parse(options.inf,options.type)

GC=[]
mult=[]
for n,read in enumerate(inf):
    if n%1000==0:
        print >> sys.stderr, 'processed {0}k reads\r'.format(n/1000),
    if options.genome is not None:
        if read.reference_name not in genome.references:
            break
        chunks=np.concatenate([[read.reference_start],read.reference_start+np.where(np.diff(read.positions) > 1)[0],[read.reference_end]])
        seq=''.join([genome.fetch(reference=read.reference_name,start=chunks[i],end=chunks[i+1]) for i in range(len(chunks)-1)]).upper()
    else:
        seq=str(read.seq).upper()
    GC.append((seq.count('G')+seq.count('C'))/float(len(seq)))
    if options.regress_multiplicity:
        mult.append(np.log10(int(read.id.split('-')[-1])))
    if options.nmax is not None and n >= options.nmax:
        break

print >> sys.stderr, 'processed {0} reads'.format(n)

GC=np.array(GC)
mult=np.array(mult)

print '# {0} ({1} reads)'.format(options.inf,len(GC))
print '# GC: mu={0:.3f} sigma={1:.3f}'.format(100*np.mean(GC),100*np.std(GC))

counts,bins=np.histogram(GC,bins=np.linspace(0,1,options.nbins+1))

if options.regress_multiplicity:
    ii=np.searchsorted(bins,GC)-1
    avgmult=[np.mean(mult[ii==i]) for i in range(options.nbins)]
    stdmult=[np.std(mult[ii==i]) for i in range(options.nbins)]
    print '# regression of log10 mult on GC: slope={0:.3f}, intercept={1:.3f}, rval={2:.3f}, pval={3:.2g}'.format(*scipy.stats.linregress(GC,mult))
    print 'GC\tfrac_reads\tavg_log10mult\tstd_log10mult'
    print '\n'.join('{0}\t{1:.6f}\t{2:.6f}\t{3:.6f}'.format(b,c/float(len(GC)),m,s) for c,b,m,s in zip(counts,bins,avgmult,stdmult))
else:
    print 'GC\tfrac_reads'
    print '\n'.join('{0}\t{1:.6f}'.format(b,c/float(len(GC))) for c,b in zip(counts,bins))

