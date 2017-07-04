import os
import sys
import numpy as np
import pandas as pd
import subprocess
import StringIO
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-i','--inf',dest='motif_definitions',help="""file with motif_definitions for logo generation""")
parser.add_option('-o','--outdir',dest='outdir',help="""output directory""")

options,args=parser.parse_args()

if not os.path.isdir(options.outdir):
	print >> sys.stderr, 'creating output directory '+options.outdir
	os.system('mkdir '+options.outdir)

print >> sys.stderr, 'creating sequence logos using motif definitions from '+options.motif_definitions
for line in open(options.motif_definitions):
	motif=line.split()[0]
	kmers=line.split()[1].split(',')
	print >> sys.stderr, motif
	seqs='\n'.join('>seq_{0}.format(n)\n{1}'.format(n,s) for n,s in enumerate(kmers))
	mout,merr=subprocess.Popen(['muscle'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate(seqs)
	outfile=os.path.join(options.outdir,'{0}.png'.format(motif))
	wout,werr=subprocess.Popen(['weblogo','-F','png_print','-A','dna','-X','no','-Y','no','-P','','-c','classic','--errorbars','no','-o',outfile],\
							   stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate(mout)
