import os
import sys
import numpy as np
import pandas as pd
import scipy.stats
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-i','--indir',dest='indir',help="""directory with deseq2 output files (from run_deseq2_for_mirca.R)""")
parser.add_option('-c','--indir_control',dest='indir_control',help="""directory with deseq2 output files for control run""")
parser.add_option('-a','--alpha',dest='alpha',default=0.05,type=float,help="""significance cutoff [0.05]""")
parser.add_option('-o','--outf',dest='outf',help="""output file [stdout]""")

options,args=parser.parse_args()

results_files=dict((f.split("deseq2_results_")[1].split(".csv")[0],f) for f in os.listdir(options.indir) if f.startswith('deseq2_results') and f.endswith('.csv'))
print >> sys.stderr, 'reading {0} results files from {1}'.format(len(results_files),options.indir)
deseq2_results=pd.concat([pd.read_csv(os.path.join(options.indir,f),index_col=0,header=0) for f in results_files.values()],axis=1,keys=results_files.keys())

if options.indir_control is not None:

	control_files=dict((f.split("deseq2_results_")[1].split(".csv")[0],f) for f in os.listdir(options.indir_control) if f.startswith('deseq2_results') and f.endswith('.csv'))
	print >> sys.stderr, 'reading {0} control files from {1}'.format(len(control_files),options.indir_control)
	deseq2_control=pd.concat([pd.read_csv(os.path.join(options.indir_control,f),index_col=0,header=0) for f in control_files.values()],axis=1,keys=control_files.keys())

	print >> sys.stderr, 'estimating p-value cutoff using control files'

	pvals=deseq2_results.xs('pvalue',axis=1,level=1).unstack().dropna().sort_values()
	pvals_c=deseq2_control.xs('pvalue',axis=1,level=1).unstack().dropna().sort_values()

	n=len(pvals)
	nc=len(pvals_c)

	ratio=np.array([float(i*nc)/float(n*np.searchsorted(pvals_c.values,p)) if p > pvals_c.iloc[0] else np.inf for i,p in enumerate(pvals.values)])
	if np.any(ratio > (1./options.alpha+1)):
		nsig = np.sum(ratio > (1./options.alpha+1))
		pc=pvals.iloc[nsig]
	else:
		nsig = 0
		pc=0

	print >> sys.stderr, "found {0} differential binding events at {1:.0f}% FDR".format(nsig,100*options.alpha)

else:

	print >> sys.stderr, 'estimating p-value cutoff using adjusted p-values from deseq2'

	pvals=deseq2_results.xs('pvalue',axis=1,level=1)
	qvals=deseq2_results.xs('padj',axis=1,level=1)

	nsig = (qvals < options.alpha).max().max()
	pc = pvals[qvals < options.alpha].max().max()

	print >> sys.stderr, "found {0} differential binding events at {1:.0f}% FDR".format(nsig,100*options.alpha)

events=deseq2_results.stack(level=0,dropna=True)
sig_events=events[events['pvalue'] < pc]
events_by_motif=sig_events.groupby(level=1)
motif_stats=events_by_motif.apply(lambda df: pd.Series(dict(ndown=(df['log2FoldChange'] < 0).sum(),\
															nup=(df['log2FoldChange'] > 0).sum(),\
															lfc_mean=df['log2FoldChange'].mean(),\
															lfc_sem=df['log2FoldChange'].sem(),\
															lfc_pval=scipy.stats.ttest_1samp(df['log2FoldChange'],0)[1])))

if options.outf is not None:
	print >> sys.stderr, 'writing to '+options.outf
	motif_stats.to_csv(options.outf)
else:
	print >> sys.stderr, 'writing to stdout'
	motif_stats.to_csv(sys.stdout)
	
