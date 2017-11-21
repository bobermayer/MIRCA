import os
import sys
import re
import numpy as np
import pandas as pd
import bisect
import gzip
import pysam
import pyBigWig
from string import maketrans
from collections import defaultdict,Counter
from optparse import OptionParser

def RC (s):

	""" get reverse complement of sequence """

	rc_tab=maketrans('ACGTUNacgtun','TGCAANtgcaan')

	return s.translate(rc_tab)[::-1]

def merge_intervals (intervals):

	""" interval merging function from here: http://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals """

	sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
	merged = []

	for higher in sorted_by_lower_bound:
		if not merged:
			merged.append(higher)
		else:
			lower = merged[-1]
			# test for intersection between lower and higher:
			# we know via sorting that lower[0] <= higher[0]
			if higher[0] <= lower[1]:
				upper_bound = max(lower[1], higher[1])
				merged[-1] = (lower[0], upper_bound)  # replace by merged interval
			else:
				merged.append(higher)

	return merged

def get_regions_from_bed(bed_file,flank_len=100):

	""" parse 12-column bed file and extract exons, CDS, UTRs and intron flanks """

	for line in bed_file:

		ls=line.split()
		chrom,tstart,tend,name,_,strand,cstart,cend,_,nexons,exon_size,exon_start=ls
		tstart,tend,cstart,cend,nexons=map(int,[tstart,tend,cstart,cend,nexons])
		exon_start=[tstart+int(x) for x in exon_start.strip(',').split(',')]
		exon_end=[exon_start[k]+int(x) for k,x in enumerate(exon_size.strip(',').split(','))]
		length=tend-tstart
		fe=bisect.bisect(exon_start,cstart)-1
		le=bisect.bisect(exon_start,cend)-1

		utr5_exons=[(exon_start[k],exon_end[k]) for k in range(fe)]+[(exon_start[fe],cstart)]
		utr3_exons=[(cend,exon_end[le])]+[(exon_start[k],exon_end[k]) for k in range(le+1,nexons)]
		if le > fe:
			cds_exons=[(cstart,exon_end[fe])]+[(exon_start[k],exon_end[k]) for k in range(fe+1,le)]+[(exon_start[le],cend)]
		else:
			cds_exons=[(cstart,cend)]
		tx_exons=[(exon_start[k],exon_end[k]) for k in range(nexons)]
		intron_up_flanks=[]
		intron_down_flanks=[]
		for k in range(nexons-1):
			intron_half_point=exon_end[k]+(exon_start[k+1]-exon_end[k])/2
			intron_up_flanks.append((exon_end[k],min(exon_end[k]+flank_len,intron_half_point)))
			intron_down_flanks.append((max(intron_half_point,exon_start[k+1]-flank_len),exon_start[k+1]))

		if strand=='-':
			utr3_exons,utr5_exons=utr5_exons,utr3_exons
			intron_up_flanks,intron_down_flanks=intron_down_flanks,intron_up_flanks

		yield (name,chrom,strand,{'utr3': utr3_exons,\
								  'utr5': utr5_exons,\
								  'cds': cds_exons,\
								  'tx': tx_exons,\
								  'intron_up': intron_up_flanks,\
								  'intron_down': intron_down_flanks})

def get_regions_from_gtf(gtf_file,flank_len=100):

	""" parse gencode-style gtf file and extract exons, CDS, UTRs and intron flanks """

	def get_exons (gene, lines, flank_len):

		""" extract regions from a bunch of lines belonging to one gene """

		all_exons=[]
		cds_exons=[]
		utr_exons=[]
		chroms=set([])
		strands=set([])

		for ls in lines:

			info=dict((x.split()[0].strip(),x.split()[1].strip().strip('"')) for x in ls[8].strip(';').split(";"))
			chrom=ls[0]
			strand=ls[6]
			
			chroms.add(chrom)
			strands.add(strand)

			# extract CDS, UTR and exon lines for this gene
			if ls[2]=='CDS':
				cds_exons.append((int(ls[3])-1,int(ls[4])))
			elif ls[2]=='UTR':
				utr_exons.append((int(ls[3])-1,int(ls[4])))
			# ignore exons of non-protein-coding isoforms of protein-coding genes
			if ls[2]=='exon' and (info['gene_type']!='protein_coding' or info['transcript_type']=='protein_coding'):
				all_exons.append((int(ls[3])-1,int(ls[4])))

		if len(chroms)!=1 or len(strands)!=1:
			raise Exception("more than one chrom or strand in lines")

		if len(cds_exons)==0:
			return (gene, chroms.pop(), strands.pop(), [])

		min_CDS=min(start for start,_ in cds_exons)
		max_CDS=max(end for _,end in cds_exons)

		utr3_exons=[]
		utr5_exons=[]
		for start,end in utr_exons:
			if end <= min_CDS:
				utr5_exons.append((start,end))
			elif start >= max_CDS:
				utr3_exons.append((start,end))

		merged_exons=merge_intervals(all_exons)

		if len(merged_exons) > 0:

			exon_start,exon_end=zip(*merged_exons)
			nexons=len(exon_start)
			intron_up_flanks=[]
			intron_down_flanks=[]
			for k in range(nexons-1):
				intron_half_point=exon_end[k]+(exon_start[k+1]-exon_end[k])/2
				intron_up_flanks.append((exon_end[k],min(exon_end[k]+flank_len,intron_half_point)))
				intron_down_flanks.append((max(intron_half_point,exon_start[k+1]-flank_len),exon_start[k+1]))

		else:

			intron_up_flanks=[]
			intron_down_flanks=[]

		if strand=='-':
			utr5_exons,utr3_exons=utr3_exons,utr5_exons
			intron_up_flanks,intron_down_flanks=intron_down_flanks,intron_up_flanks

		return (gene, chroms.pop(), strands.pop(),{'utr3':merge_intervals(utr3_exons),\
												   'utr5':merge_intervals(utr5_exons),\
												   'cds':merge_intervals(cds_exons),\
												   'tx':merged_exons,\
												   'intron_up':intron_up_flanks,\
												   'intron_down':intron_down_flanks})

	gene_lines=defaultdict(list)
	for line in gtf_file:

		if line.startswith('#'):
			continue

		ls=line.strip().split("\t")
		info=dict((x.split()[0].strip(),x.split()[1].strip().strip('"')) for x in ls[8].strip(';').split(";"))
		name=info['gene_id']
		gene_lines[name].append(ls)

	for gene,lines in gene_lines.iteritems():
		yield get_exons(gene,lines,flank_len)

if __name__ == '__main__':

	parser=OptionParser()
	parser.add_option('-i','--inf',dest='inf',help="either: 12-column bed file with transcript definitions; or GTF file")
	parser.add_option('-b','--bigWig',dest='bigWig',help="bigWig file with conservation values")
	parser.add_option('-K','--K',dest='K',default=7,help="k-mers to analyze [7]",type=int)
	parser.add_option('-E','--E',dest='E',default=0,help="extend k-mer regions by E nucleotides [0]",type=int)
	parser.add_option('-f','--flank_len',dest='flank_len',default=100,type=int,help="length of intron flanks [100]")
	parser.add_option('-r','--region',dest='region',default='utr3',help="which region to use (utr5/cds/utr3/tx/intron_up/intron_down) [utr3]")
	parser.add_option('-m','--motif_file',dest='motif_file',help="file with motif name and comma-separated list of kmers for each motif")
	parser.add_option('-g','--genome',dest='genome',help="genome file (fasta; must be indexed)")
	parser.add_option('-o','--outfile',dest='outf',help="output file [stdout]")

	options,args=parser.parse_args()

	K=options.K
	E=options.E

	if options.outf is None:
		outf=sys.stdout
	else:
		outf=open(options.outf,'w')

	if options.region not in ['utr5','cds','utr3','tx','intron_up','intron_down']:
		raise Exception("unknown region selected!")

	print >> sys.stderr, 'conserved motifs in {0} from {1}'.format(options.region,options.inf)
	outf.write('# conserved motifs in {0} from {1}\n'.format(options.region,options.inf))

	use_motifs=False
	if options.motif_file is not None:
		print >> sys.stderr, 'using RBP motifs from '+options.motif_file
		outf.write('# using RBP motifs from {0}\n'.format(options.motif_file))
		motifs=pd.DataFrame([(k,line.split()[0]) for line in open(options.motif_file) for k in line.split()[1].split(',')],columns=['kmer','motif']).set_index('kmer')
		use_motifs=True
	print >> sys.stderr, 'K={0}, E={1}, flank_len={2}'.format(K,E,options.flank_len)
	outf.write('# K={0}, E={1}, flank_len={2}\n'.format(K,E,options.flank_len))

	genome=pysam.FastaFile(options.genome)

	if options.bigWig is not None:
		print >> sys.stderr, 'using bigWig file',options.bigWig
		bw=pyBigWig.open(options.bigWig)
		outf.write('# bigWig file: {0}\n'.format(options.bigWig))

	outf.write('gene\t{0}\tcount\tcons\n'.format('motif' if use_motifs else 'kmer'))

	if '.gtf' in options.inf:
		get_regions = get_regions_from_gtf
	elif '.bed' in options.inf:
		get_regions = get_regions_from_bed
	else:
		raise Exception("cannot recognize format of {0}!".format(options.inf))

	if options.inf.endswith('.gz'):
		inf=gzip.open(options.inf,'rb')
	else:
		inf=open(options.inf)

	nskipped=0

	for name,chrom,strand,regions in get_regions(inf, options.flank_len):

		if len(regions)==0:
			continue

		region_exons=regions[options.region]
		if len(region_exons)==0:
			continue

		cons=[]
		seqs=[]

		for start,end in region_exons:
			if end <= start:
				continue
			# get conservation
			cons.append(np.nan_to_num(bw.values(chrom,start,end)))
			# get sequence for the entire region 
			seqs.append(genome.fetch(reference=chrom,start=start,end=end).upper())

		if len(cons)==0 or len(seqs)==0:
			continue

		# combine exons for utr3/utr5/cds/tx
		if 'intron' not in options.region:
			cons=[np.concatenate(cons,axis=0)]
			seqs=[''.join(seqs)]

		# reverse orientation if on minus strand
		if strand=='-':
			cons=map(lambda x: x[::-1], cons[::-1])
			seqs=map(RC, seqs[::-1])

		exon_length=map(len,seqs)
		nexons=len(seqs)

		if sum(map(len,cons))!=sum(exon_length):
			raise Exception("lengths don't match!")

		# don't consider this gene or transcript if there is no sequence
		if sum(exon_length) < K:
			nskipped+=1
			continue

		# get count and conservation for each kmer in each region
		kmer_cons=[]
		for i in range(nexons):
			for k in range(exon_length[i]-K+1):
				# use mean read count over stretch of length K+2*E
				kmer_cons.append([seqs[i][k:k+K],1,cons[i][max(k-E,0):min(k+K+E,exon_length[i])].mean()])
			kmer_cons.append(['tot',sum(exon_length),cons[i].sum()])

		# sum kmer counts and average conservation over all occurrences of the kmer
		tot_kmer_cons=pd.DataFrame(kmer_cons, columns=['kmer','count','cons']).groupby('kmer').sum()

		# print motif occurrences and conservation
		if use_motifs:
			# sum over all kmers for this motif
			tot_motif_cons=tot_kmer_cons.join(motifs).groupby('motif').sum()
			tot_motif_cons.loc['tot']=tot_kmer_cons.loc['tot']
			for motif,vals in tot_motif_cons.iterrows():
				outf.write("{0}\t{1}\t{2:.0f}\t{3:.4f}\n".format(name,motif,vals['count'],vals['cons']/vals['count']))
		else:
			for kmer,vals in tot_kmer_cons.iterrows():
				outf.write("{0}\t{1}\t{2:.0f}\t{3:.4f}\n".format(name,kmer,vals['count'],vals['cons']/vals['count']))

	print >> sys.stderr, 'done ({0} skipped)'.format(nskipped)

	if outf is not sys.stdout:
		outf.close()
