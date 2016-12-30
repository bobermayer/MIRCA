# MIRCA

MIRCA (Motif-Informed Read Coverage Analysis) is a collection of scripts to find differentially bound RBP motif occurrences in protein occupancy profiling data (see [Baltz et al. (2012)](http://dx.doi.org/10.1016/j.molcel.2012.05.021)). In contrast to [POPPI](http://dx.doi.org/10.1186/gb-2014-15-1-r15), which looks for differentially bound *positions*, here we aggregate counts of reads or T->C conversions over all occurrences of predefined motifs (e.g., taken from [CISBP-RNA](http://cisbp-rna.ccbr.utoronto.ca)) and then check for differences between conditions in entire transcripts or transcript regions using counts over the entire region as covariate for expression normalization with [DESeq2](http://dx.doi.org/10.1186/s13059-014-0550-8).

## Prerequisites
MIRCA uses Python (including numpy, scipy, pandas, pysam, and Biopython), R (including optparse, DESeq2, and data.table), and samtools. Required inputs are bam files with mapped reads, transcript annotation (gencode-style GTF or 12-column bed), genome fasta, and a file with motif definitions (see the attached file Mouse_RNAcompete_kmers_condensed.txt as an example; it was prepared by selecting for each motif the Kmers with Z-scores at least 5 sigma over the mean and then combining motifs with at least 90% overlap between their Kmers).

## Description
The first script (get_mirca_read_counts.py) extracts regions (all exons, UTRs, CDS or intron flanks) from the GTF or bed file and uses samtools mpileup to count reads (or T->C conversion events) over all words of length K (Kmers) in the region for each bam file supplied. Word occurrences can be extended by E nucleotides to capture neighboring conversion events. Counts for Kmers belonging to the same motif are added up. The second script (run_deseq2_for_mirca.R) then uses DESeq2 to assess differential occupancy (relative to read or conversion counts over the entire region) between conditions. More complex designs can be tested by editing a few lines in the script. This yields statistics for each motif in each region, which are collected using the script collect_mirca_results.py. Significant events are filtered either using DESeq2 adjusted p-value (which are very conservative), or more sensitively by comparing to control results with permuted labels (similar to POPPI).

## Usage

### 1. get read counts
```
python get_mirca_read_counts.py \
       -i input.bed/input.gtf(.gz) \
       -B bam_condition1_rep1,bam_condition1_rep2,bam_condition2_rep1,bam_condition2_rep2 \
       -K 7 \
       -E 0 \
       -r utr5/cds/utr3/tx/intron_up/intron_down \
       -m rbp_motifs.txt \
       -g genome.fa > mirca_counts.out
```
Conversion events are counted when using ``-T``. Names for files can be given with ``-n "condition1_1,condition1_2,condition2_1,condition2_2"``, otherwise bam file names are used. Additional options can be explored using ``python get_mirca_read_counts.py -h``.

### 2. run DESeq2
``` Rscript run_deseq2_for_mirca.R -i mirca_counts.out -o mirca_deseq2_results -c condition1,condition1,condition2,condition2 ```
Here, ``-o`` is an output directory (will be created if not exists) with deseq2 results for each motif, and ``-c`` specifies the conditions corresponding to the columns of mirca_counts.out. By default, DESeq2 will test differential occupancy between condition2 and condition1, but more complex designs can be tested by appropriately editing lines 57-61 in the script.

A control run with permuted labels can be performed like so:
``` Rscript run_deseq2_for_mirca.R -i mirca_counts.out -o mirca_deseq2_control -c condition1,condition2,condition1,condition2 ```


### 3. collect results
``` python collect_mirca_results.py -i mirca_deseq2_results -o mirca_results.tsv ```
This script collects the DESeq2 output in mirca_deseq2_results and selects significant events using adjusted p-values. An empirical FDR for selecting significant events (specified with ``-a``, default: 0.05) can be calculated by specifying the output directory of a control run using ``-c mirca_deseq2_control``. Events per motif are combined and statistics collected: number of genes with up-regulated binding, number of genes with down-regulated binding, mean and s.e.m of occupancy log2 fold change across genes, a p-value from a t-test of these log2 fold changes against 0, and lists of up- and down-regulated genes. 
