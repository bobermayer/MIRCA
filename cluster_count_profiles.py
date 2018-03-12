import os
import sys
import numpy as np
import pandas as pd
import scipy.spatial
import scipy.cluster
import fastcluster
from sklearn.decomposition import PCA
from optparse import OptionParser

def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

parser=OptionParser()
parser.add_option('-i','--counts_in',dest='counts_in',\
				  help="""input file with read counts (output of get_mirca_counts.py)""")
parser.add_option('-o','--counts_out',dest='counts_out',default="clustered_counts.out",\
				  help="""output file with clustered read counts [clustered_counts.out]""")
parser.add_option('','--method',dest='method',default='ward',\
				  help="""clustering method (one of single, complete, average, weighted, ward, centroid, median) [ward]""")
parser.add_option('','--metric',dest='metric',default='correlation',\
				  help="""distance metric (used by fastcluster) [correlation]""")
parser.add_option('','--t',dest='t',\
				  help="""cutoff for scipy.cluster.hierarchy.fcluster [0.7*max(Z[:,2])]""")
parser.add_option('','--criterion',dest='criterion',default='distance',\
				  help="""cutoff criterion for scipy.cluster.hierarchy.fcluster [distance]""")
parser.add_option('','--quantile',dest='quantile',default=0.1,type=float,\
				  help="""quantile cutoff on coverage to select genes used for clustering [0.1]""")
parser.add_option('','--max_motifs_per_cluster',dest='max_motifs_per_cluster',default=500,type=int,\
				  help="""max num of motifs per [500]""")
parser.add_option('','--use_pca',dest='use_pca',default=False,action='store_true',\
				  help="""use dimensionality reduction (PCA) before clustering (not recommended)""")
parser.add_option('','--use_variable',dest='use_variable',default=False,action='store_true',\
				  help="""use only variable motif-gene combinations before clustering""")
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

motifs=np.setdiff1d(counts.index.get_level_values(1).unique(),['tot'])
conditions=counts.columns[1:]
	
# normalize counts by "tot" counts over regions
tot_counts=counts[conditions].xs('tot',level=1,axis=0).astype(float)

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
	for n,mot in enumerate(motifs):
		print >> sys.stderr, 'normalizing counts for {0} of {1} motifs\r'.format(n+1,len(motifs)),
		normalized_counts[n]=counts_here.loc[mot].divide(tot_counts,axis=0).stack().reindex(ii,fill_value=0).values
	print >> sys.stderr, ''

# use only motifs for which there is data
take=np.sum(normalized_counts,axis=1) > 0
motifs=motifs[take]
normalized_counts=normalized_counts[take]

if options.use_variable:
	# use only variable gene-motif combinations for clustering (based on mean-CV dependence)
	means=normalized_counts.mean(axis=0)
	CV=normalized_counts.std(axis=0)/means

	ok=(means > np.percentile(means,5)) & (means < np.percentile(means,95))
	ii=np.argsort(means[ok])
	CV_mean=np.mean(rolling_window(CV[ok][ii],100),-1)
	CV_std=np.std(rolling_window(CV[ok][ii],100),-1)
	zscores=(CV[ok][ii][50:-49]-CV_mean)/CV_std
	take=np.arange(normalized_counts.shape[1])[ok][ii][50:-49][zscores > .5]
	take.sort()
	normalized_counts=normalized_counts[:,take]

if options.use_pca:
	pca=PCA(n_components=50,copy=False,whiten=False).fit(normalized_counts.T)
	components=pd.DataFrame(pca.components_,index=np.arange(50),columns=motifs)
	ncomponents=np.min(np.where(pca.explained_variance_ratio_[:-1]/pca.explained_variance_ratio_[1:] < 1.05)[0])
	print >> sys.stderr, 'reducing data to {0} PCs'.format(ncomponents)
	data=components.ix[:ncomponents].T
else:
	data=normalized_counts

# cluster count profiles over genes
Z0=fastcluster.linkage(data,metric=options.metric,method=options.method)
threshold=0.7*max(Z0[:,2]) if options.t is None else float(options.t)
clusters0=scipy.cluster.hierarchy.fcluster(Z0,threshold,criterion=options.criterion)

# retain only clusters with less than options.max_matifs_per_cluster of motifs
good_clusters=[n for n in np.unique(clusters0) if np.sum(clusters0==n) < options.max_motifs_per_cluster]
take=np.in1d(clusters0,good_clusters)
motifs=motifs[take]
Z=fastcluster.linkage(data[take],metric=options.metric,method=options.method)
threshold=0.7*max(Z[:,2]) if options.t is None else float(options.t)
clusters=scipy.cluster.hierarchy.fcluster(Z,threshold,criterion=options.criterion)

# get motifs for each cluster
cluster_motifs=dict(('cluster_{0}'.format(c),list(motifs[clusters==c])) for c in np.unique(clusters))
nclusters=len(cluster_motifs)
print >> sys.stderr, 'obtained {0} clusters after dropping {1}'.format(nclusters,len(set(clusters0)-set(good_clusters)))

print >> sys.stderr, 'combining counts over clusters'
# get combined dataframe of counts: sum up counts over all motifs belonging to one cluster and add total counts
clustered_counts=pd.concat([counts[counts.index.get_level_values(1).isin(clust)].sum(axis=0,level=0) \
							for clust in cluster_motifs.values()] + \
                               [counts.xs('tot',level=1,axis=0).astype(float)],\
						   axis=0,keys=cluster_motifs.keys()+['tot']).swaplevel(0,1,axis=0).sort_index(axis=0).fillna(0).astype(int)[counts.columns]
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
	from collections import Counter
	import subprocess

	matplotlib.rcParams['lines.linewidth'] = 0.5

	plt.ion()

	if options.motif_definitions is not None:
		print >> sys.stderr, 'reading motif definitions from '+options.motif_definitions
		motif_kmers=pd.read_csv(options.motif_definitions,header=None,index_col=0,sep='\t').squeeze()
		cluster_seqs=dict((n,'\n'.join('>cluster_{0}.seq_{1}\n{2}'.format(n+1,k+1,m) for k,m in enumerate(set(','.join(motif_kmers[motifs[clusters==(n+1)]].dropna()).split(','))))) for n in range(nclusters))
	else:
		cluster_seqs=dict((n,'\n'.join('>cluster_{0}.seq_{1}\n{2}'.format(n+1,k+1,m) for k,m in enumerate(motifs[clusters==(n+1)]))) for n in range(nclusters))

	print >> sys.stderr, 'getting cluster consensus sequences'
	cluster_labels=[]
	for n in range(nclusters):
		mout,merr=subprocess.Popen(['muscle'],\
								   stdin=subprocess.PIPE,\
								   stdout=subprocess.PIPE,\
								   stderr=subprocess.PIPE).communicate(cluster_seqs[n])
		maf=np.array([list(line) for line in mout.split('\n') if not line.startswith('>') and len(line) > 0])
		consensus=''
		for k in range(maf.shape[1]):
			col=Counter(maf[:,k])
			c=col.keys()[np.argmax(col.values())]
			if c=='-':
				consensus+='-'
			elif col[c] >= .75*len(maf):
				consensus+=c.upper()
			elif col[c] >= .25*len(maf):
				consensus+=c.lower()
			else:
				consensus+='.'
		cluster_labels.append('cluster_{0}: {1}'.format(n+1,consensus))

	fig=plt.figure(figsize=(7,12))
	fig.clf()

	ax=fig.add_axes([.02,.05,.4,.9])
	zz=scipy.cluster.hierarchy.dendrogram(Z,orientation='right',labels=motifs,ax=ax)
	#zz=scipy.cluster.hierarchy.dendrogram(Z,orientation='left',no_labels=True,ax=ax)
	for t in ax.get_yticklabels():
		t.set_fontsize(6)

	ax.set_xticks([])
	#ax.set_xlim(sorted(ax.get_xlim())[::-1])
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)

	cluster_edges=np.concatenate([[0],np.where(np.diff(clusters[zz['leaves']])!=0)[0]+1,[len(clusters)]])
	cluster_ids=[np.unique(clusters[zz['leaves'][cluster_edges[n]:cluster_edges[n+1]]])[0] for n in range(nclusters)]
	for n in range(nclusters):
		k=cluster_ids[n]-1
		pos_x=-2
		pos_y=10*(n+.5)*len(motifs)/float(nclusters)
		pos_y2=5.*(cluster_edges[n]+cluster_edges[n+1])
		ax.text(pos_x,pos_y,
				cluster_labels[k],size=6,\
				va='center',family='monospace',ha='left',color='grcmyk'[k%6],\
				bbox=dict(facecolor='none',edgecolor='grcmyk'[k%6]),clip_on=False)
		ax.vlines(pos_x+.5,10*cluster_edges[n]+5,10*cluster_edges[n+1]-5,color='grcmyk'[k%6],lw=.5,clip_on=False)
		ax.arrow(pos_x+.5,pos_y2,-.4,pos_y-pos_y2,
				 color='grcmyk'[k%6],lw=.5,head_width=.0,head_length=.0,clip_on=False)

	print >> sys.stderr, 'saving figure '+options.out_fig
	fig.savefig(options.out_fig)

