library('optparse')
library(DESeq2)
library('data.table')

option_list <- list(
  make_option(c("-i","--input"),type="character",help="input file (from get_mirca_read_counts.py)",metavar="file"),
  make_option(c("-o","--outdir"),type='character',default="mirca_out",help="output dir [%default]",metavar="dir"),
  make_option(c("-c","--conditions"),type="character",help="comma-separated list of conditions (corresponding to the columns in input file)",metavar="list"),
  make_option("--alpha",type="double",default=0.01,help="significance cutoff alpha [%default]",metavar="double")
)

opt_parser <- OptionParser(option_list = option_list)
options <- parse_args(opt_parser)

dir.create(file.path(options$outdir), showWarnings = TRUE)

if (is.null(options$input)) {
  print_help(opt_parser)
  stop("please supply input file",call.=FALSE)
}

if (is.null(options$conditions)) {
  print_help(opt_parser)
  stop("please supply comma-separated list of conditions",call.=FALSE)
}

# ADD COVARIATES FOR MORE COMPLEX DESIGN IF AVAILABLE
conditions <- strsplit(options$conditions,',',fixed=TRUE)[[1]]

print(paste("reading from",options$input))
counts.all <- read.csv(options$input,sep='\t',comment='#',header=T)

sampleName <- colnames(counts.all)[3:length(colnames(counts.all))]
if (length(conditions) != length(sampleName)) {
  stop("number of conditions (",length(conditions),") doesn't match number of columns (",length(sampleName),") in ",options$input,call.=FALSE)
}
# take only named conditions, drop empty ones or those called "_"
take <- (conditions!='_') & (conditions!='')
colData <- data.frame(condition=conditions[take],row.names=sampleName[take])
columns <- c('gene',sampleName[take])
print(paste('using columns',paste(columns[2:length(columns)],collapse=','),'for conditions',paste(conditions[take],collapse=',')))

all.motifs <- setdiff(unique(counts.all$motif),'tot')

counts.tot <- data.table(na.omit(counts.all[counts.all$motif=='tot',columns]))
counts.tot.df <- as.data.frame.matrix(counts.tot[,2:length(counts.tot),with=FALSE],row.names=as.vector(counts.tot$gene))
dds.tot <- DESeqDataSetFromMatrix(counts.tot.df,colData,~ condition)
sf.tot <- estimateSizeFactorsForMatrix(counts(dds.tot))

for (mot in all.motifs) {
  
    print (paste('run deseq2 for',mot))
    
    counts.motif <- data.table(na.omit(counts.all[counts.all$motif==mot,columns]))
    counts.here <- merge(counts.tot,counts.motif,by='gene',all.y=TRUE,suffixes=c('_tot',paste('_',mot,sep='')))
    counts.here.df <- as.data.frame.matrix(counts.here[,2:length(counts.here),with=FALSE],row.names=as.vector(counts.here$gene))
    sampleName <- colnames(counts.here)[2:length(colnames(counts.here))]
    motifs <- c(rep("tot",sum(take)),rep(mot,sum(take)))
    # EDIT COLDATA AND FULL/REDUCED DESIGN IF MORE COVARIATES ARE AVAILABLE
    colData <- data.frame(condition=c(conditions[take],conditions[take]),motif=motifs,row.names=sampleName)
    dds <- DESeqDataSetFromMatrix(counts.here.df, colData, ~ condition + motif + condition:motif)
    dds$motif <- factor(dds$motif, levels=c('tot',mot))
    sizeFactors(dds) <- c(sf.tot,sf.tot)
    tryCatch({
      dds <- DESeq(dds, test='LRT', reduced=~condition+motif)
      res <- results(dds,alpha=options$alpha)
      write.csv(as.data.frame(res),file=file.path(options$outdir,paste('deseq2_results_',mot,'.csv',sep='')))
    }, error=function(err) {
      print(paste('skipping',mot,'due to error:',err))
    })
}  

