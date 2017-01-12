# MIRCA

MIRCA (Motif-Informed Read Coverage Analysis) is a collection of scripts to find differentially bound RBP motif occurrences in protein occupancy profiling data (see [Baltz et al. (2012)](http://dx.doi.org/10.1016/j.molcel.2012.05.021)). In contrast to [POPPI](http://dx.doi.org/10.1186/gb-2014-15-1-r15), which looks for differentially bound *positions*, here we aggregate counts of reads or T->C conversions over all occurrences of predefined motifs (e.g., taken from [CISBP-RNA](http://cisbp-rna.ccbr.utoronto.ca)) in entire transcripts or transcript regions. We then check for differences between conditions using [DESeq2](http://dx.doi.org/10.1186/s13059-014-0550-8), using counts over the entire region for expression normalization.

## Prerequisites
MIRCA uses Python (including numpy, scipy, pandas, matplotlib, and pysam) and R (including optparse, DESeq2, and data.table). Required inputs are bam files with mapped reads, transcript annotation (gencode-style GTF or 12-column bed), indexed genome fasta, and a file with motif definitions (see the attached [file](Mouse_RNAcompete_kmers_condensed.txt) as an example; it was prepared using mouse motifs and Zscores from [CISBP-RNA](http://cisbp-rna.ccbr.utoronto.ca) by selecting for each motif the Kmers with Z-scores at least 5 sigma over the mean and then combining motifs with at least 90% overlap between their Kmers). If this motif file is not given, read/conversion counts are aggregated over all Kmers in the transcriptome separately. Clustering count profiles for each motif or Kmer requires fastcluster, Biopython, subprocess and clustalw2.

## Description
The first script ``get_mirca_read_counts.py`` extracts regions (all exons, UTRs, CDS or intron flanks) for each gene (GTF file) or transcript (bed file) and uses samtools mpileup to count reads or conversions over all words of length K (Kmers) in the region for each bam file supplied. Word occurrences can be extended by E nucleotides to capture neighboring conversion events, which are often observed nearby but not within miRNA seed matches, for instance. If a file with motif definitions is given, counts for Kmers belonging to the same motif are added up. Alternatively, count profiles for different Kmers (or motifs) can be clustered using ``cluster_count_profiles.py``, which reduces the number of Kmers/motifs and therefore reduces running time and increases statistical power of downstream analysis. However, this is a somewhat experimental feature requiring manual tuning by the user (as most clustering routines do). The second script ``run_deseq2_for_mirca.R`` then uses DESeq2 to assess differential occupancy (relative to read or conversion counts over the entire region) between conditions. More complex designs can be tested by editing a few lines in the script. This yields statistics for each motif in each region, which are collected using another script ``collect_mirca_results.py``. Significant events are filtered either using DESeq2 adjusted p-value (which are very conservative), or more sensitively by comparing to control results with permuted labels (similar to POPPI).

## Usage

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

### 2. (optional) cluster count profiles for different motifs or Kmers
```
python cluster_count_profiles.py \
	-i mirca_counts.out \
	-o mirca_counts_clustered.out
```
Count profiles over genes and conditions for clusters of motifs (or Kmers) defined via hierarchical clustering are added up. The clustering metric and method can be specified with ``--metric`` and ``--method``, respectively, as well as other arguments to ``scipy.cluster.hierarchy.fcluster``. If an argument is provided to ``--fig``, this script also plots the associated dendrogram and does a very simple calculation of cluster consensus sequences by aligning with ``clustalw2`` (if motifs are clustered, motif definitions used for ``get_mirca_read_counts.py`` need to be provided with ``--motif_definitions``). The top 40% of genes (by total counts) can be selecting ``--quantile 0.6``.

**Note**: This script can take a lot of memory and computing time (linear in the number of genes times number of conditions, and quadratic in the number of motifs/Kmers), and will likely require some manual parameter tuning by the user to provide meaningful cluster selection.

### 3. run DESeq2
``` 
Rscript run_deseq2_for_mirca.R \
	-i mirca_counts.out \
	-o mirca_deseq2_results \
	-c condition1,condition1,condition2,condition2 
```
Here, ``-o`` is an output directory (will be created if it doesn't exists) with deseq2 results for each motif, and ``-c`` specifies the conditions corresponding to the columns of ``mirca_counts.out``. Zero-length conditions or conditions called ``_`` will be ignored, allowing to select specific columns in the input file. DESeq2 is used with a likelihood ratio test (LRT) for the significance of an interaction between the condition and a ``motif`` factor, which selects read/conversion counts over a specific motif or the entire region, respectively (i.e., changes in the ratio of read counts over motif vs. the entire regions are tested for differences between conditions; by default, differential occupancy is tested between the alphabetically last and the first condition in this list, but different or more complex designs can be tested by appropriately editing lines 60-61 in the script).

A control run with permuted labels can be performed like this:
`` Rscript run_deseq2_for_mirca.R -i mirca_counts.out -o mirca_deseq2_control -c condition1,condition2,condition1,condition2 ``


### 4. collect results
``` 
python collect_mirca_results.py -i mirca_deseq2_results -o mirca_results.tsv 
```

This script collects the DESeq2 output in ``mirca_deseq2_results`` and selects significant events. An empirical FDR for selecting significant events (specified with ``-a``, default: 0.05) can be calculated by specifying the output directory of a control run using ``-c mirca_deseq2_control``, otherwise DESeq2 adjusted p-values are used. Events per motif are combined and statistics collected: number of genes with up-regulated binding, number of genes with down-regulated binding, fraction of differentially bound targets (diff. bound genes relative to all genes with nonzero coverage at this motif), mean and s.e.m. of occupancy log2 fold change across differentially bound genes, a p-value from a t-test of these log2 fold changes against 0, mean and s.e.m. of occupancy log 2 fold change across all genes, and lists of up- and down-regulated genes. 

If ``-f output_plot.pdf -t plot_title`` are added, this script also produces a summary plot.
