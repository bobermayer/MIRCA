# MIRCA

MIRCA (Motif-Informed Read Coverage Analysis) is a collection of scripts to find differentially bound RBP motif occurrences in protein occupancy profiling data (see [Baltz et al. (2012)](http://dx.doi.org/10.1016/j.molcel.2012.05.021)). In contrast to [POPPI](http://dx.doi.org/10.1186/gb-2014-15-1-r15), which looks for differentially bound *positions*, here we aggregate counts of reads or T->C conversions over all occurrences of predefined motifs (e.g., taken from [CISBP-RNA](http://cisbp-rna.ccbr.utoronto.ca)) in entire transcripts or transcript regions. We then check for differences between conditions using [DESeq2](http://dx.doi.org/10.1186/s13059-014-0550-8), using counts over the entire region for expression normalization.

## Prerequisites
MIRCA uses Python (including numpy, scipy, pandas, matplotlib, and pysam) and R (including optparse, DESeq2, and data.table). Required inputs are bam files with mapped reads, transcript annotation (gencode-style GTF or 12-column bed), indexed genome fasta, and a file with motif definitions (see the attached [file](Mouse_RNAcompete_kmers_condensed.txt) as an example; it was prepared using mouse motifs and Zscores from [CISBP-RNA](http://cisbp-rna.ccbr.utoronto.ca) by selecting for each motif the Kmers with Z-scores at least 5 sigma over the mean and then combining motifs with at least 90% overlap between their Kmers). Clustering count profiles for each motif or Kmer requires fastcluster, Biopython, subprocess and muscle. Logo creation is also supported if [weblogo](https://pypi.python.org/pypi/weblogo) is in the path.

## Description
The first script ``get_mirca_read_counts.py`` extracts regions (all exons, UTRs, CDS or intron flanks) for each gene (GTF file) or transcript (bed file) and uses samtools mpileup to count reads or conversions over all words of length K (Kmers) in the region for each bam file supplied. Word occurrences can be extended by E nucleotides to capture neighboring conversion events, which are often observed nearby but not within miRNA seed matches, for instance. If a file with motif definitions is given, counts for Kmers belonging to the same motif are added up. Alternatively, count profiles for different Kmers (or motifs) can be clustered using ``cluster_count_profiles.py``, which reduces the number of Kmers/motifs and therefore reduces running time and increases statistical power of downstream analysis. However, this is a somewhat experimental feature requiring manual tuning by the user (as most clustering routines do). Additionally, the script simply counts motif occurrences in the selected sequence regions (may be useful for enrichment calculations). The second script ``run_deseq2_for_mirca.R`` then uses DESeq2 to assess differential occupancy (relative to read or conversion counts over the entire region) between conditions. More complex designs can be tested by editing a few lines in the script. This yields statistics for each motif in each region, which are collected using another script ``collect_mirca_results.py``. Significant events are filtered either using adjusted p-values from DESeq2, or more permissively by comparing to control results with permuted labels (similar to POPPI).

## Usage

### 0. (optional) balance GC content of input reads
differences in GC content between libraries can lead to spurious associations with GC- or AT-rich motifs. The following two scripts allow to characterize the GC distribution in a bam file and then sample from these files to match a target distribution.  

``` 
python get_read_GC_content.py -i input.bam > input_GC_stats.txt
python get_GC_balanced_reads.py -i input.bam -o input_balanced.bam \
	--target_GC target_GC_stats.txt --input_GC input_GC_stats.txt > input_balanced_GC_stats.txt
```
the target GC distribution can be obtained, e.g., by averaging over different input files:
```
paste input_1_GC_stats.txt input_2_GC_stats.txt ... | grep -v "^#" |  \
	awk '{s=n=0; for (i=2; i<=NF; i+=2) {s+=$i;n++} printf "%s\t%.6f\n",$1,(n ? s/n : $0)}' > target_GC_stats.txt
```
fastq (or fastq.gz) files can also be used in both scripts by specifying ``-t fastq``.

### 1. get read counts
```
python get_mirca_read_counts.py \
       -i input.bed/input.gtf(.gz) \
       -B condition1_rep1.bam,condition1_rep2.bam,condition2_rep1.bam,condition2_rep2.bam \
       -K 7 \
       -r utr5/cds/utr3/tx/intron_up/intron_down \
       -m rbp_motifs.txt \
       -g genome.fa \
       -o mirca_counts.out
```
Conversion events instead of reads are counted when using ``-T``. More readable names for files can be given with ``-n "condition1_1,condition1_2,condition2_1,condition2_2"``, otherwise bam file names are used. Additional options can be explored using ``python get_mirca_read_counts.py -h``, e.g., the intron flank length can be adjusted using ``-f`` (default: 100nt), and Kmer windows can be extended by ``n`` nucleotides using ``-E n``. 

### 2. (optional) cluster count profiles 
```
python cluster_count_profiles.py \
	-i mirca_counts.out \
	-o mirca_counts_clustered.out
```
Count profiles over genes and conditions for clusters of motifs (or Kmers) defined via hierarchical clustering are added up. Cluster members are written within in comments in the header section of the output file, which would then be used in step 3 instead of ``mirca_counts.out``. The clustering metric and method can be specified with ``--metric`` and ``--method``, respectively, as well as other arguments to ``scipy.cluster.hierarchy.fcluster``. If an argument is provided to ``--fig``, this script also plots the associated dendrogram and does a very simple calculation of cluster consensus sequences by aligning with ``muscle`` (if motifs containing Kmers are clustered, motif definitions used for ``get_mirca_read_counts.py`` need to be provided with ``--motif_definitions``). The top 40% of genes (by total counts) can be selecting ``--quantile 0.6``.

**Note**: This script can take a lot of memory and computing time (linear in the number of genes times number of conditions, and quadratic in the number of motifs/Kmers), and will likely require some manual parameter tuning by the user to give meaningful clusters.

### 3. run DESeq2
``` 
Rscript run_deseq2_for_mirca.R \
	-i mirca_counts.out \
	-o mirca_results.csv.gz \
	-c condition1,condition1,condition2,condition2 
```
Here, ``-c`` specifies the conditions corresponding to the columns of ``mirca_counts.out``. Zero-length conditions or conditions called ``_`` will be ignored, allowing to select specific columns in the input file. By default, all conditions will be tested against a reference condition (alphabetically first one, or one given using ``--reference``). DESeq2 is used with a likelihood ratio test (LRT) for the significance of an interaction between the condition and a ``motif`` factor, which selects read/conversion counts over a specific motif or the entire region, respectively (i.e., changes in the ratio of read counts over motif vs. the entire regions are tested for differences between conditions. Different or more complex designs can be tested by appropriately editing the script at the indicated lines. Numerical or categorical covariates for each condition can be supplied with ``--num_covariate`` or ``--cat_covariate``.

A control run with permuted labels can be performed like this:
`` Rscript run_deseq2_for_mirca.R -i mirca_counts.out -o mirca_control.csv.gz -c condition1,condition2,condition1,condition2 ``

### 4. collect results
``` 
python collect_mirca_results.py -i mirca_results.csv.gz -s mirca_results_summary.tsv 
```

This script uses the DESeq2 output in ``mirca_results.csv.gz`` and determines globally significant events. An empirical FDR for selecting significant events (specified with ``-a``, default: 0.05) can be calculated by specifying the output of a control run using ``-c mirca_control.csv.zg``, otherwise p-values for all genes and all motifs are globally adjusted using the BH method. 

If ``-s`` is given, then events per motif are combined and summary statistics collected: number of genes with up-regulated binding, number of genes with down-regulated binding, fraction of differentially bound targets (diff. bound genes relative to all genes with nonzero coverage at this motif), mean and s.e.m. of occupancy log2 fold change across differentially bound genes, a p-value from a t-test of these log2 fold changes against 0, mean and s.e.m. of occupancy log 2 fold change across all genes, and lists of up- and down-regulated genes. 

If ``-f output_plot.pdf -t plot_title`` are added, this script also produces a summary plot. If additionally ``--motif_definitions motif_definitions.txt`` is given, the script also adds logos for each motif or cluster to the plot (motif definitions from clustered counts can be extracted using ``grep "^# cluster_" clustered_counts.out | tr -d "# " | tr ":" "\t" > motif_definitions.txt``).

