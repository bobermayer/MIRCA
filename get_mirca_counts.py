import os
import sys
import numpy as np
import pandas as pd
import bisect
import gzip
import pysam
import subprocess
from Bio import SeqIO
from string import maketrans
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

	def get_exons (lines, flank_len):

		""" extract regions from a bunch of lines belonging to one gene """

		all_exons=[]
		cds_exons=[]
		utr_exons=[]
		names=set([])
		chroms=set([])
		strands=set([])

		for line in lines:

			ls=line.strip().split('\t')

			info=dict((x.split()[0].strip(),x.split()[1].strip().strip('"')) for x in ls[8].strip(';').split(";"))
			name=info['gene_id']
			chrom=ls[0]
			strand=ls[6]
			
			names.add(name)
			chroms.add(chrom)
			strands.add(strand)

			if ls[2]=='gene':
				raise Exception("this shouldn't happen!")
			elif ls[2] in ['CDS','start_codon','stop_codon']:
				cds_exons.append((int(ls[3])-1,int(ls[4])))
			elif ls[2]=='UTR':
				utr_exons.append((int(ls[3])-1,int(ls[4])))
			# ignore exons of non-protein-coding isoforms of protein-coding genes
			if ls[2]=='exon' and (info['gene_type']!='protein_coding' or info['transcript_type']=='protein_coding'):
				all_exons.append((int(ls[3])-1,int(ls[4])))

		if len(names)!=1 or len(chroms)!=1 or len(strands)!=1:
			raise Exception("more than one gene or chrom or strand in lines")

		if len(cds_exons)==0:
			return (names.pop(),chroms.pop(),strands.pop(),[])

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

		return (names.pop(),chroms.pop(),strands.pop(),{'utr3':merge_intervals(utr3_exons),\
														'utr5':merge_intervals(utr5_exons),\
														'cds':merge_intervals(cds_exons),\
														'tx':merged_exons,\
														'intron_up':intron_up_flanks,\
														'intron_down':intron_down_flanks})

	lines=[]
	for line in gtf_file:

		if line.startswith('#'):
			continue

		ls=line.strip().split("\t")

		if ls[2]=='gene':
			if len(lines) > 0:
				yield get_exons(lines,flank_len)
				lines=[]
		else:
			lines.append(line)

	if len(lines) > 0:
		yield get_exons(lines,flank_len)
		lines=[]

if __name__ == '__main__':

	parser=OptionParser()
	parser.add_option('-i','--inf',dest='inf',help="either: 12-column bed file with transcript definitions; or GTF file")
	parser.add_option('-B','--bam',dest='bam',help="comma-separated list of BAM files with mapped reads, should have indices")
	parser.add_option('-K','--K',dest='K',default=7,help="k-mers to analyze (default: 7)",type=int)
	parser.add_option('-E','--E',dest='E',default=5,help="extend k-mer regions by E nucleotides (default: 5)",type=int)
	parser.add_option('-f','--flank_len',dest='flank_len',default=100,type=int,help="length of intron flanks (default: 100)")
	parser.add_option('-r','--region',dest='region',default='utr3',help="which region to use (utr5/cds/utr3/tx/intron_up/intron_down) (default: utr3)")
	parser.add_option('-n','--names',dest='names',default=None,help="header names for bam files (comma-separated)")
	parser.add_option('-m','--motif_file',dest='motif_file',help="file with motif name and comma-separated list of kmers for each motif")
	parser.add_option('-g','--genome',dest='genome',help="genome file (fasta; must be indexed)")
	parser.add_option('-T','--use_TC',dest='use_TC',action='store_true',help="use number of TC conversions instead of read counts")

	options,args=parser.parse_args()

	K=options.K
	E=options.E

	if options.region not in ['utr5','cds','utr3','tx','intron_up','intron_down']:
		raise Exception("unknown region selected!")

	use_motifs=False
	if options.motif_file is not None:
		print >> sys.stderr, 'using RBP motifs from '+options.motif_file
		motifs=pd.DataFrame([(k,line.split()[0]) for line in open(options.motif_file) for k in line.split()[1].split(',')],columns=['kmer','motif']).set_index('kmer')
		use_motifs=True

	print >> sys.stderr, 'using bam files',options.bam
	bam_files=[bam.strip() for bam in options.bam.split(',')]
	nmapped=np.array([pysam.Samfile(bam,'rb').mapped for bam in bam_files])
	nB=len(bam_files)

	print >> sys.stderr, 'counts in {0} from '.format(options.region)+options.inf

	if options.names is not None:
		names=dict((n,x.strip()) for n,x in enumerate(options.names.split(',')))
		if len(names)!=nB:
			raise Exception("number of header names doesn't match number of bam files")
	else:
		names=dict(zip(range(nB),range(1,nB+1)))

	sys.stdout.write('# {0} counts in {1} from {2}'.format('TC' if options.use_TC else 'read', options.region,options.inf)+'\n# bam files:\n')
	for n in range(nB):
		sys.stdout.write('#  {0}: {1} ({2} reads)\n'.format(names[n],options.bam.split(',')[n],nmapped[n]))

	if options.motif_file is not None:
		sys.stdout.write('# using RBP motifs from {0}\n'.format(options.motif_file))
	sys.stdout.write('# K={0}, E={1}, flank_len={2}'.format(K,E,options.flank_len))

	sys.stdout.write('gene\t{0}'.format('motif' if use_motifs else 'kmer'))
	for n in range(nB):
		sys.stdout.write('\tcov_{0}'.format(names[n]))
	sys.stdout.write('\n')

	if 'gtf' in options.inf:
		get_regions = get_regions_from_gtf
	else:
		get_regions = get_regions_from_bed

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

		covs=[]
		seqs=[]

		for start,end in region_exons:
			# get read coverage using samtools mpileup
			proc=subprocess.Popen(['samtools','mpileup']+bam_files+['-f',options.genome,'-r',chrom+':'+str(start+1)+'-'+str(end)],stdout=subprocess.PIPE)
			cov=np.zeros((end-start,nB),dtype=np.int)
			while True:
				line=proc.stdout.readline()
				if line=='':
					break
				ls=line.strip('\n').split('\t')
				pos=int(ls[1])-1
				ref=ls[2]
				if pos < start or pos > end-1:
					raise Exception("unknown position in mpileup")
				if options.use_TC:
					if strand=='+' and ref in 'Tt':
						cov[pos-start]=np.array([ls[3*(n+1)+1].count('C')+ls[3*(n+1)+1].count('c') for n in range(nB)])
					elif strand=='-' and ref in 'Aa':
						cov[pos-start]=np.array([ls[3*(n+1)+1].count('G')+ls[3*(n+1)+1].count('g') for n in range(nB)])
				else:
					cov[pos-start]=np.array([int(ls[3*(n+1)]) for n in range(nB)])
			covs.append(cov)
			# get sequence for the entire region (even where there is no coverage)
			proc=subprocess.Popen(['samtools','faidx']+[options.genome,chrom+':'+str(start+1)+'-'+str(end)],stdout=subprocess.PIPE)
			seq=str(SeqIO.read(proc.stdout,'fasta').seq).upper()
			seqs.append(seq)

		# combine exons for utr3/utr5/cds/tx
		if 'intron' not in options.region:
			covs=[np.concatenate(covs,axis=0)]
			seqs=[''.join(seqs)]

		# reverse orientation if on minus strand
		if strand=='-':
			covs=map(lambda x: x[::-1], covs[::-1])
			seqs=map(RC, seqs[::-1])

		exon_length=map(len,seqs)
		nexons=len(seqs)

		if sum(map(len,covs))!=sum(exon_length):
			raise Exception("lengths don't match!")

		# don't consider this gene or transcript if there is no coverage
		if sum(exon_length) <= 0 or sum(map(lambda x: x.sum(),covs))==0:
			nskipped+=1
			continue

		# get coverage for each kmer in each region
		kmer_cov=[]
		for i in range(nexons):
			for k in range(exon_length[i]-K):
				if np.all(covs[i][k:k+K].sum(axis=0)==0):
					continue
				if options.use_TC:
					# use actual number of T->C conversions in that region
					kmer_cov.append([seqs[i][k:k+K]]+[covs[i][max(k-E,0):min(k+K+E,exon_length[i]),n].sum() for n in range(nB)])
				else:
					# use mean read count over stretch of length K+2*E
					kmer_cov.append([seqs[i][k:k+K]]+[covs[i][max(k-E,0):min(k+K+E,exon_length[i]),n].mean() for n in range(nB)])
			kmer_cov.append(['tot']+list(covs[i].sum(axis=0)))

		if len(kmer_cov) > 0:
			# sum kmer coverage over all occurrences of the kmer
			tot_kmer_cov=pd.DataFrame(kmer_cov, columns=['kmer']+map(lambda x: 'cov_{0}'.format(names[x]),range(nB))).groupby('kmer').sum()
			if use_motifs:
				# sum over all kmers for this motif
				tot_motif_cov=tot_kmer_cov.join(motifs).groupby('motif').sum()
				tot_motif_cov.loc['tot']=tot_kmer_cov.loc['tot']
				# print counts as integer
				for motif,vals in tot_motif_cov.iterrows():
					sys.stdout.write("{0}\t{1}\t".format(name,motif)+'\t'.join('{0:.0f}'.format(v) for v in vals)+'\n')
			else:
				# print counts as integer
				for kmer,vals in tot_kmer_cov.iterrows():
					sys.stdout.write("{0}\t{1}\t".format(name,kmer)+'\t'.join('{0:.0f}'.format(v) for v in vals)+'\n')

	print >> sys.stderr, 'done ({0} skipped)'.format(nskipped)