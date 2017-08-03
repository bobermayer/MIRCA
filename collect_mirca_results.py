import os
import sys
import gzip
import numpy as np
import pandas as pd
import scipy.stats
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def collect_stats (df, alpha):

	""" helper function to aggregate binding events for each motif """

	# all significant events
	sig=df['qvalue'] < alpha
	# all upregulated events
	up=df['log2FoldChange'] > 0
	# all downregulated events
	down=df['log2FoldChange'] < 0
	# all bound events (with coverage > 0)
	bound=df['baseMean'] > 0

	return pd.Series([(up & sig).sum(),\
					  (down & sig).sum(),\
					  sig.sum()/float(bound.sum()),\
					  df[sig]['log2FoldChange'].mean(),\
					  df[sig]['log2FoldChange'].sem(),\
					  scipy.stats.ttest_1samp(df[sig]['log2FoldChange'],0)[1],\
					  df['log2FoldChange'].mean(),\
					  df['log2FoldChange'].sem(),\
					  ','.join(zip(*df[sig & up].index.tolist())[0]) if (sig & up).sum() > 0 else '',\
					  ','.join(zip(*df[sig & down].index.tolist())[0]) if (sig & down).sum() > 0 else ''],\
					 index=['nup_sig','ndown_sig','frac_sig','lfc_mean_sig','lfc_sem_sig','lfc_pval_sig','lfc_mean_all','lfc_sem_all','genes_up','genes_down'])

def p_adjust_bh(p):

	""" Benjamini-Hochberg p-value correction for multiple hypothesis testing
		(taken and modified for nan values from here: http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python) """

	p = np.asfarray(p)
	ok = np.isfinite(p)
	by_descend = p[ok].argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(p[ok])) / np.arange(len(p[ok]), 0, -1)
	q = np.zeros_like(p)
	q[ok] = np.minimum(1, np.minimum.accumulate(steps * p[ok][by_descend]))[by_orig]
	q[~ok] = np.nan
	return q

def get_empirical_qval (pvals,pvals_c):

	""" get empirical q-values by comparing p-values against control p-values with permutated labels """

	p=np.asfarray(pvals)
	pc=np.asfarray(pvals_c.unstack().dropna().sort_values())

	ok=np.isfinite(p)
	by_ascend = p[ok].argsort()
	by_orig = by_ascend.argsort()
	q = np.zeros_like(p)
	q[ok] = 2./(1.+np.arange(len(p[ok])).astype(float)/np.searchsorted(pc,p[ok][by_ascend]))[by_orig]
	q[~ok] = np.nan

	return q

parser=OptionParser()
parser.add_option('-i','--inf',dest='inf',help="""deseq2 output file (from run_deseq2_for_mirca.R)""")
parser.add_option('-c','--inf_control',dest='inf_control',help="""deseq2 output file for control run""")
parser.add_option('-a','--alpha',dest='alpha',default=0.05,type=float,help="""significance cutoff [0.05]""")
parser.add_option('-s','--summary',dest='summary',help="""summary output [None]""")
parser.add_option('-f','--fig',dest='fig',help="""figure output [None]""")
parser.add_option('-t',dest='title',default='MIRCA stats',help="""title for figure""")
parser.add_option('','--motif_images',dest='motif_images',help="""folder with logos for motifs (created by make_logos.py)""")

options,args=parser.parse_args()

print >> sys.stderr, 'reading deseq2 results files from {0}'.format(options.inf)
deseq2_results=pd.read_csv(options.inf,index_col=0,header=0)
deseq2_results.columns=pd.MultiIndex.from_tuples([c.split('.') for c in deseq2_results.columns])

if options.inf_control is not None:

	print >> sys.stderr, 'reading deseq2 control files from {0}'.format(options.inf_control)
	deseq2_control=pd.read_csv(options.inf_control,index_col=0,header=0)
	deseq2_control.columns=pd.MultiIndex.from_tuples([c.split('.') for c in deseq2_control.columns])

	print >> sys.stderr, 'estimating empirical q-values using control files'

	pvals=deseq2_results.xs('pvalue',axis=1,level=2)
	pvals_c=deseq2_control.xs('pvalue',axis=1,level=2)

	qvals = pd.DataFrame(get_empirical_qval (pvals, pvals_c),index=pvals.index,columns=pvals.columns)

else:

	print >> sys.stderr, 'estimating q-values using global BH correction'

	pvals=deseq2_results.xs('pvalue',axis=1,level=2)
	qvals = pd.DataFrame(p_adjust_bh(pvals),index=pvals.index,columns=pvals.columns)

print >> sys.stderr, "found {0} differential binding events at {1:.0f}% FDR".format((qvals < options.alpha).sum().sum(),100*options.alpha)

motifs=deseq2_results.columns.get_level_values(0).unique()
conditions=deseq2_results.columns.get_level_values(1).unique()

# add these q-values to the table
for mot in motifs:
	for cond in conditions:
		deseq2_results[mot,cond,'qvalue']=qvals[mot,cond]
deseq2_results=deseq2_results.sort_index(axis=1)

# combine events by motif and condition and collect stats
combined_results=deseq2_results.stack(level=[0,1],dropna=True).groupby(level=[1,2]).apply(lambda df: collect_stats (df,options.alpha)).unstack(level=1).swaplevel(0,1,axis=1)

if options.summary is not None:

	print >> sys.stderr, 'writing summary to '+options.summary
	combined_results.to_csv(options.summary,sep='\t')
	
if options.fig is not None:

	totdiff=combined_results.xs('nup_sig',axis=1,level=1).sum(axis=1)+\
			 combined_results.xs('ndown_sig',axis=1,level=1).sum(axis=1)

	order=totdiff.sort_values().index

	figheight=.15*len(order)+1
	figwidth=2.+4*len(conditions)

	left=2./figwidth
	dist=4/figwidth
	bottom=.4/figheight
	height=1-.35/figheight-bottom
	width=.8*4/figwidth/3.
	wspace=.05/len(conditions)

	if options.motif_images is not None:

		from matplotlib import image as mpimg
		print >> sys.stderr, 'reading motif images'
		motif_images=dict((mot,mpimg.imread(os.path.join(options.motif_images,'{0}.png'.format(mot)))) for mot in order)

	fig=plt.figure(1,figsize=(figwidth,figheight))
	fig.clf()

	for n,cond in enumerate(conditions):

		ax=fig.add_axes([left+n*dist,bottom,width,height])
		ax.barh(np.arange(len(order)),combined_results.ix[order,(cond,'nup_sig')],color='r',label='up',lw=0)
		ax.barh(np.arange(len(order)),combined_results.ix[order,(cond,'ndown_sig')],left=combined_results.ix[order,(cond,'nup_sig')],color='b',label='down',lw=0)
		ax.set_yticks(np.arange(len(order))+.4)
		ax.set_yticklabels([o if len(o) < 15 else o[:12]+'...' for o in order] if n==0 else [],size=8,ha='right',va='center')
		ax.set_xlabel('diff. bound genes',size=8)
		ax.set_ylim([0,len(order)-.2])
		ax.set_xlim([0,1.1*max(map(abs,ax.get_xlim()))])
		ax.locator_params(axis='x',nbins=5,integer=True)
		ax.tick_params(axis='x', which='major', labelsize=8)
		ax.grid(axis='x')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().tick_left()
		#ax.get_xaxis().set_major_locator(MaxNLocator(integer=True))

		if n==0:
			leg=ax.legend(loc=4,prop={'size':8},title='occupancy')
			leg.get_title().set_fontsize(8)
			leg.get_frame().set_ec('k')

		ax=fig.add_axes([left+n*dist+width+wspace,bottom,width,height])
		ax.barh(np.arange(len(order)),100*combined_results.ix[order,(cond,'frac_sig')],color='g',lw=0)
		ax.set_yticks(np.arange(len(order))+.4)
		ax.set_yticklabels([])
		ax.set_xlabel('% targets',size=8)
		ax.set_ylim([0,len(order)-.2])
		ax.set_xlim([0,1.1*max(ax.get_xlim())])
		ax.locator_params(axis='x',nbins=3)
		ax.tick_params(axis='x', which='major', labelsize=8)
		ax.grid(axis='x')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().tick_left()
		if len(conditions) > 1:
			ax.set_title(cond,size=8)

		ax=fig.add_axes([left+n*dist+2*(width+wspace),bottom,width,height])
		ax.barh(np.arange(len(order)),combined_results.ix[order,(cond,'lfc_mean_sig')],xerr=combined_results.ix[order,(cond,'lfc_sem_sig')],color='m',height=.8,label='diff',ecolor='m',error_kw={'lw':.5},lw=0)
		ax.barh(np.arange(len(order))+.1,combined_results.ix[order,(cond,'lfc_mean_all')],color='DarkGray',height=.6,label='all',lw=0)
		ax.set_yticks(np.arange(len(order))+.4)
		ax.set_yticklabels([])
		ax.set_xlabel('occupancy log2 FC',size=8)
		ax.set_ylim([0,len(order)-.2])
		ax.locator_params(axis='x',nbins=5)
		ax.tick_params(axis='x', which='major', labelsize=8)
		ax.grid(axis='x')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().tick_left()

		if n==0:
			leg=ax.legend(loc=4,prop={'size':8},title='genes')
			leg.get_title().set_fontsize(8)
			leg.get_frame().set_ec('k')

	for n,motif in enumerate(order):
		if motif in motif_images:
			ax=fig.add_axes([.01/len(conditions),bottom+n*height/len(order),.5*left,height/len(order)])
			ax.imshow(motif_images[motif])
			ax.set_axis_off()
		
	fig.suptitle(options.title+' ({0:.0f}% FDR)'.format(100*options.alpha),size=8,y=.7+.3*(bottom+height),va='center')

	print >> sys.stderr, 'saving figure to '+options.fig
	fig.savefig(options.fig,dpi=600)
