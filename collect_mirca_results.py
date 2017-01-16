import os
import sys
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
parser.add_option('-i','--indir',dest='indir',help="""directory with deseq2 output files (from run_deseq2_for_mirca.R)""")
parser.add_option('-c','--indir_control',dest='indir_control',help="""directory with deseq2 output files for control run""")
parser.add_option('-a','--alpha',dest='alpha',default=0.05,type=float,help="""significance cutoff [0.05]""")
parser.add_option('-o','--outf',dest='outf',help="""output file [stdout]""")
parser.add_option('-s','--summary',dest='summary',help="""summary output [None]""")
parser.add_option('-f','--fig',dest='fig',help="""figure output [None]""")
parser.add_option('-t',dest='title',default='MIRCA stats',help="""title for figure""")

options,args=parser.parse_args()

results_files=dict((f.split("deseq2_results_")[1].split(".csv")[0],f) for f in os.listdir(options.indir) if f.startswith('deseq2_results') and f.endswith('.csv'))
print >> sys.stderr, 'reading {0} results files from {1}'.format(len(results_files),options.indir)
deseq2_results=pd.concat([pd.read_csv(os.path.join(options.indir,f),index_col=0,header=0) for f in results_files.values()],axis=1,keys=results_files.keys())

if options.indir_control is not None:

	control_files=dict((f.split("deseq2_results_")[1].split(".csv")[0],f) for f in os.listdir(options.indir_control) if f.startswith('deseq2_results') and f.endswith('.csv'))
	print >> sys.stderr, 'reading {0} control files from {1}'.format(len(control_files),options.indir_control)
	deseq2_control=pd.concat([pd.read_csv(os.path.join(options.indir_control,f),index_col=0,header=0) for f in control_files.values()],axis=1,keys=control_files.keys())

	print >> sys.stderr, 'estimating empirical q-values using control files'

	pvals=deseq2_results.xs('pvalue',axis=1,level=1)
	pvals_c=deseq2_control.xs('pvalue',axis=1,level=1)

	qvals = pd.DataFrame(get_empirical_fdr (pvals, pvals_c),index=pvals.index,columns=pvals.columns)

else:

	print >> sys.stderr, 'estimating q-values using global BH correction'

	pvals=deseq2_results.xs('pvalue',axis=1,level=1)
	qvals = pd.DataFrame(p_adjust_bh(pvals),index=pvals.index,columns=pvals.columns)

print >> sys.stderr, "found {0} differential binding events at {1:.0f}% FDR".format((qvals < options.alpha).sum().sum(),100*options.alpha)

# add these q-values to the table
for m in deseq2_results.columns.get_level_values(0).unique():
	deseq2_results[m,'qvalue']=qvals[m]
deseq2_results=deseq2_results.sort_index(axis=1)

# combine events by motif and collect stats
events_by_motif=deseq2_results.stack(level=0,dropna=True).groupby(level=1)
motif_stats=events_by_motif.apply(lambda df: collect_stats (df,options.alpha))

if options.summary is not None:

	print >> sys.stderr, 'writing summary to '+options.summary
	motif_stats.to_csv(options.summary,sep='\t')
	
if options.fig is not None:

	order=(motif_stats['nup_sig']+motif_stats['ndown_sig']).sort_values().index

	fig=plt.figure(1,figsize=(8,10))
	fig.clf()

	ax=fig.add_axes([.37,.08,.2,.89])
	ax.barh(np.arange(len(order)),motif_stats.ix[order,'nup_sig'],color='r',label='up',lw=0)
	ax.barh(np.arange(len(order)),motif_stats.ix[order,'ndown_sig'],left=motif_stats.ix[order,'nup_sig'],color='b',label='down',lw=0)
	ax.set_yticks(np.arange(len(order))+.4)
	ax.set_yticklabels(order,size=6,ha='right',va='center')
	ax.set_xlabel('diff. bound genes',size=8)
	ax.set_ylim([-.4,len(order)+.2])
	ax.set_xlim([0,1.1*max(map(abs,ax.get_xlim()))])
	ax.locator_params(axis='x',nbins=5)
	ax.tick_params(axis='x', which='major', labelsize=8)
	ax.grid(axis='x')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	leg=ax.legend(loc=3,prop={'size':8},ncol=2,bbox_to_anchor=(0,-.09),title='occupancy')
	leg.get_title().set_fontsize(8)
	leg.get_frame().set_ec('k')

	ax=fig.add_axes([.59,.08,.18,.89])
	ax.barh(np.arange(len(order)),100*motif_stats.ix[order,'frac_sig'],color='g',lw=0)
	ax.set_yticks(np.arange(len(order))+.4)
	ax.set_yticklabels([])
	ax.set_xlabel('% targets',size=8)
	ax.set_ylim([-.4,len(order)+.2])
	ax.set_xlim([0,1.1*max(ax.get_xlim())])
	ax.locator_params(axis='x',nbins=3)
	ax.tick_params(axis='x', which='major', labelsize=8)
	ax.grid(axis='x')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()

	ax=fig.add_axes([.8,.08,.18,.89])
	ax.barh(np.arange(len(order)),motif_stats.ix[order,'lfc_mean_sig'],xerr=motif_stats.ix[order,'lfc_sem_sig'],color='m',height=.8,label='diff',ecolor='m',error_kw={'lw':.5},lw=0)
	ax.barh(np.arange(len(order))+.1,motif_stats.ix[order,'lfc_mean_all'],color='DarkGray',height=.6,label='all',lw=0)
	ax.set_yticks(np.arange(len(order))+.4)
	ax.set_yticklabels([])
	ax.set_xlabel('occupancy log2 FC',size=8)
	ax.set_ylim([-.4,len(order)+.2])
	ax.locator_params(axis='x',nbins=5)
	ax.tick_params(axis='x', which='major', labelsize=8)
	ax.grid(axis='x')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	leg=ax.legend(loc=3,prop={'size':8},ncol=2,bbox_to_anchor=(0,-.09),title='genes')
	leg.get_title().set_fontsize(8)
	leg.get_frame().set_ec('k')

	fig.suptitle(options.title+' ({0:.0f}% FDR)'.format(100*options.alpha),size=8,y=.99)

	print >> sys.stderr, 'saving figure to '+options.fig
	fig.savefig(options.fig)

# flatten hierarchical index and write output
deseq2_results.columns=['.'.join(c) for c in deseq2_results.columns.tolist()]
print >> sys.stderr, 'writing full output to '+(options.outf if options.outf is not sys.stdout else 'stdout')
deseq2_results.to_csv(options.outf,sep='\t')
