import os
import sys
import numpy as np
import pandas as pd
import scipy.spatial
import scipy.cluster
import fastcluster
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-i','--counts_in',dest='counts_in',\
				  help="""input file with read counts (output of get_mirca_counts.py)""")
parser.add_option('-o','--counts_out',dest='counts_out',default="clustered_counts.out",\
				  help="""output file with clustered read counts [clustered_counts.out]""")
parser.add_option('','--method',dest='method',default='ward',\
				  help="""clustering method (one of single, complete, average, weighted, ward, centroid, median) [ward]""")
parser.add_option('','--metric',dest='metric',default='correlation',\
				  help="""distance metric (as in scipy.distance.pdist) [correlation]""")
parser.add_option('','--t',dest='t',\
				  help="""cutoff for scipy.cluster.hierarchy.fcluster [0.7*max(Z[:,2])]""")
parser.add_option('','--criterion',dest='criterion',default='distance',\
				  help="""cutoff criterion for scipy.cluster.hierarchy.fcluster [distance]""")
parser.add_option('','--quantile',dest='quantile',default=0.1,type=float,\
				  help="""quantile cutoff on coverage to select genes used for clustering [0.1]""")
parser.add_option('','--motif_definitions',dest='motif_definitions',
				  help="""file with motif_definitions used in get_mirca_read_counts.py (otherwise, motif names are used to get consensus)""")
parser.add_option('','--fig',dest='out_fig',\
				  help="""output figure with dendrogram and consensus motifs""")

options,args=parser.parse_args()

# get read counts from input file
print >> sys.stderr, 'reading counts from '+options.counts_in
with open(options.counts_in) as inf:
	comments=[]
	while True:
		line=inf.readline()
		if not line.startswith('#'):
			columns=line.strip('\n').split('\t')
			break
		else:
			comments.append(line)
	counts=pd.read_csv(inf,header=None,index_col=[0,1],names=columns,sep='\t')

#################
### EDIT HERE ###
# motifs=np.random.choice(np.setdiff1d(counts.index.get_level_values(1).unique(),['tot']),1000)
#################
motifs=np.setdiff1d(counts.index.get_level_values(1).unique(),['tot'])
conditions=counts.columns
	
# normalize counts by "tot" counts over regions
tot_counts=counts.xs('tot',level=1,axis=0).astype(float)

# use only top genes for clustering
genes=tot_counts.index[tot_counts.mean(axis=1) > tot_counts.mean(axis=1).quantile(options.quantile)]
counts_here=counts[counts.index.get_level_values(0).isin(genes) &\
				   counts.index.get_level_values(1).isin(motifs)].swaplevel(0,1).sort_index(axis=0)

motifs=counts_here.index.get_level_values(0).unique()
genes=counts_here.index.get_level_values(1).unique()

# assemble normalized counts into big matrix
if counts_here.shape[0] > .1*len(motifs)*len(genes):
	print >> sys.stderr, 'normalizing counts'
	# easier way (but needs A LOT of memory when most motif x gene combinations do not appear)
	normalized_counts=counts_here.divide(tot_counts,axis=0,level=1).unstack(level=1).fillna(0).values
else:
	# slow but more memory-efficient alternative (almost equivalent except if some condition x gene combinations do not appear)
	normalized_counts=np.zeros((len(motifs),len(genes)*len(conditions)),dtype=np.float32)
	ii=pd.MultiIndex.from_product([genes,conditions])
	maxlen=max(map(len,motifs))
	for n,mot in enumerate(motifs):
		print >> sys.stderr, 'normalizing counts for {0} of {1} motifs\r'.format(n+1,len(motifs)),
		normalized_counts[n]=counts_here.loc[mot].divide(tot_counts,axis=0).stack().reindex(ii,fill_value=0).values
	print >> sys.stderr, ''

print >> sys.stderr, 'computing linkage: ',
# cluster count profiles over genes
Z=fastcluster.linkage(normalized_counts,metric=options.metric,method=options.method)
threshold=0.7*max(Z[:,2]) if options.t is None else float(options.t)
clusters=scipy.cluster.hierarchy.fcluster(Z,threshold,criterion=options.criterion)
# get motifs for each cluster
cluster_motifs=dict(('cluster_{0}'.format(c),list(motifs[clusters==c])) for c in np.unique(clusters))
nclusters=len(cluster_motifs)
print >> sys.stderr, '{0} clusters'.format(nclusters)

print >> sys.stderr, 'combining counts over clusters'
# get combined dataframe of counts: sum up counts over all motifs belonging to one cluster and add total counts
clustered_counts=pd.concat([counts[counts.index.get_level_values(1).isin(clust)].sum(axis=0,level=0) \
							for clust in cluster_motifs.values()] + [tot_counts],\
						   axis=0,keys=cluster_motifs.keys()+['tot']).swaplevel(0,1,axis=0).sort_index(axis=0).fillna(0).astype(int)
clustered_counts.index.names=['gene','motif']

print >> sys.stderr, 'writing cluster definitions and clustered counts to '+options.counts_out
with open(options.counts_out,'w') as outf:
	outf.write(''.join(comments))
	outf.write('# clustered over {0} genes, with {1} metric, {2} method, t={3:.2f} and {4} criterion\n'.format(len(genes),options.metric,options.method,threshold,options.criterion))
	for clust,clust_motifs in cluster_motifs.iteritems():
		outf.write('# {0}: {1}\n'.format(clust,','.join(clust_motifs)))
	clustered_counts.to_csv(outf,sep='\t')

if options.out_fig is not None:

	import matplotlib
	matplotlib.use('Agg')

	from matplotlib import pyplot as plt
	from Bio import SeqIO,Seq,SeqRecord,AlignIO
	from collections import Counter
	import subprocess

	matplotlib.rcParams['lines.linewidth'] = 0.5

	plt.ion()

	if options.motif_definitions is not None:
		print >> sys.stderr, 'reading motif definitions from '+options.motif_definitions
		motif_kmers=pd.read_csv(options.motif_definitions,header=None,index_col=0,sep='\t').squeeze()
		cluster_seqs=dict((n,[SeqRecord.SeqRecord(Seq.Seq(m),id='cluster_{0}.seq_{1}'.format(n+1,k+1)) for k,m in enumerate(set(','.join(motif_kmers[motifs[clusters==(n+1)]].dropna()).split(',')))]) for n in range(nclusters))
	else:
		cluster_seqs=dict((n,[SeqRecord.SeqRecord(Seq.Seq(m),id='cluster_{0}.seq_{1}'.format(n+1,k+1)) for k,m in enumerate(motifs[clusters==(n+1)])]) for n in range(nclusters))

	print >> sys.stderr, 'getting cluster consensus sequences'
	cluster_labels=[]
	for n in range(nclusters):
		SeqIO.write(cluster_seqs[n],'tmp.fa','fasta')
		tmp=subprocess.Popen(['clustalw2','-infile=tmp.fa','-align','-type=DNA','-output=PHYLIP','-outorder=INPUT','-PWGAPOPEN=10000'],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
		maf=AlignIO.read('tmp.phy','phylip')
		consensus=''
		for k in range(maf.get_alignment_length()):
			col=Counter(maf[:,k])
			c=col.keys()[np.argmax(col.values())]
			if c=='-':
				consensus+='-'
			elif col[c] >= .6*len(maf):
				consensus+=c.upper()
			elif col[c] >= .2*len(maf):
				consensus+=c.lower()
			else:
				consensus+='.'
		cluster_labels.append(consensus)

	fig=plt.figure(figsize=(7,12))
	fig.clf()

	ax=fig.add_axes([.05,.05,.4,.9])
	zz=scipy.cluster.hierarchy.dendrogram(Z,orientation='left',labels=motifs,ax=ax)

	ax.set_xticks([])
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	cluster_edges=np.concatenate([[0],np.where(np.diff(clusters[zz['leaves']])==1)[0]+1,[len(clusters)]])
	for n in range(nclusters):
		fig.text(.85,.05+.9*(cluster_edges[n]+cluster_edges[n+1])/(2.*len(clusters)),\
				 cluster_labels[n],size=6,\
				 va='center',family='monospace',ha='center',color='grcmyk'[n%6],\
				 bbox=dict(facecolor='none',edgecolor='grcmyk'[n%6]))

	print >> sys.stderr, 'saving figure '+options.out_fig
	fig.savefig(options.out_fig)

