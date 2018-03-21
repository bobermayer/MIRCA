import subprocess
from get_mirca_counts import *

if __name__ == '__main__':

        parser=OptionParser()
        parser.add_option('-i','--inf',dest='inf',help="either: 12-column bed file with transcript definitions; or GTF file")
        parser.add_option('-K','--K',dest='K',default=7,help="k-mers to analyze [7]",type=int)
        parser.add_option('-E','--E',dest='E',default=0,help="extend k-mer regions by E nucleotides [0]",type=int)
        parser.add_option('-f','--flank_len',dest='flank_len',default=100,type=int,help="length of intron flanks [100]")
        parser.add_option('-r','--region',dest='region',default='utr3',help="which region to use (utr5/cds/utr3/tx/intron_up/intron_down/intron) [utr3]")
        parser.add_option('-m','--motif_file',dest='motif_file',help="file with motif name and comma-separated list of kmers for each motif")
        parser.add_option('-g','--genome',dest='genome',help="genome file (fasta; must be indexed)")
        parser.add_option('-o','--outfile',dest='outf',help="output file [stdout]")
        parser.add_option('-t','--tmpdir',dest='tmpdir',default='.',help="tmpdir [.]")

        options,args=parser.parse_args()

        K=options.K
        E=options.E

        if options.outf is None:
                outf=sys.stdout
        else:
                outf=open(options.outf,'w')

        if options.region not in ['utr5','cds','utr3','tx','intron_up','intron_down','intron']:
                raise Exception("unknown region selected!")

        print >> sys.stderr, 'openness for motifs in {0} from {1}'.format(options.region,options.inf)
        outf.write('# openness for motifs in {0} from {1}\n'.format(options.region,options.inf))
        print >> sys.stderr, 'using temporary directory: {0}'.format(options.tmpdir)

        use_motifs=False
        if options.motif_file is not None:
                print >> sys.stderr, 'using RBP motifs from '+options.motif_file
                outf.write('# using RBP motifs from {0}\n'.format(options.motif_file))
                motifs=pd.DataFrame([(k,line.split()[0]) for line in open(options.motif_file) for k in line.split()[1].split(',')],columns=['kmer','motif']).set_index('kmer')
                use_motifs=True
        print >> sys.stderr, 'K={0}, E={1}, flank_len={2}'.format(K,E,options.flank_len)
        outf.write('# K={0}, E={1}, flank_len={2}\n'.format(K,E,options.flank_len))

        genome=pysam.FastaFile(options.genome)

        outf.write('gene\t{0}\tcount\topenness\n'.format('motif' if use_motifs else 'kmer'))

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

                openness=[]
                seqs=[]

                for start,end in region_exons:
                        if end <= start:
                                continue
                        # get sequence for the entire region 
                        seqs.append(genome.fetch(reference=chrom,start=start,end=end).upper())

                if len(seqs)==0:
                        continue

                # combine exons for utr3/utr5/cds/tx
                if 'intron' not in options.region:
                        seqs=[''.join(seqs)]

                # reverse orientation if on minus strand
                if strand=='-':
                        seqs=map(RC, seqs[::-1])

                # get openness
                for seq in seqs:
                	proc=subprocess.Popen(['RNAplfold'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=False,cwd=options.tmpdir)
                        output=proc.communicate(input=seq)[0]
                        read_probs=False
                        opness=np.ones(len(seq))
                        with open(os.path.join(options.tmpdir,"plfold_dp.ps")) as inf:
                                for line in inf:
                                        if line.startswith('%start of base pair probability data'):
                                                read_probs=True
                                        if read_probs and 'ubox' in line:
                                                i,j=map(lambda x: int(x)-1,line.split()[:2])
                                                p=float(line.split()[2])**2
                                                if i >=0 and i < len(seq):
                                                        opness[i]-=p
                                                if j >=0 and j < len(seq):
                                                        opness[j]-=p
                        openness.append(np.maximum(opness,0))                 

                exon_length=map(len,seqs)
                nexons=len(seqs)

                if sum(map(len,openness))!=sum(exon_length):
                        raise Exception("lengths don't match!")

                # don't consider this gene or transcript if there is no sequence
                if sum(exon_length) < K:
                        nskipped+=1
                        continue

                # get count and openness for each kmer in each region
                kmer_open=[]
                for i in range(nexons):
                        for k in range(exon_length[i]-K+1):
                                # use mean read count over stretch of length K+2*E
                                kmer_open.append([seqs[i][k:k+K],1,openness[i][max(k-E,0):min(k+K+E,exon_length[i])].mean()])
                        kmer_open.append(['tot',sum(exon_length),openness[i].sum()])

                # sum kmer counts and average openness over all occurrences of the kmer
                tot_kmer_open=pd.DataFrame(kmer_open, columns=['kmer','count','openness']).groupby('kmer').sum()

                # print motif occurrences and openness
                if use_motifs:
                        # sum over all kmers for this motif
                        tot_motif_open=tot_kmer_open.join(motifs).groupby('motif').sum()
                        tot_motif_open.loc['tot']=tot_kmer_open.loc['tot']
                        for motif,vals in tot_motif_open.iterrows():
                                outf.write("{0}\t{1}\t{2:.0f}\t{3:.4f}\n".format(name,motif,vals['count'],vals['openness']/vals['count']))
                else:
                        for kmer,vals in tot_kmer_open.iterrows():
                                outf.write("{0}\t{1}\t{2:.0f}\t{3:.4f}\n".format(name,kmer,vals['count'],vals['openness']/vals['count']))

        print >> sys.stderr, 'done ({0} skipped)'.format(nskipped)

        if outf is not sys.stdout:
                outf.close()
