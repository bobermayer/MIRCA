suppressPackageStartupMessages({
  require(optparse)
  require(DESeq2)
  require(data.table)
})

option_list <- list(
  make_option(c("-i","--input"),type="character",help="input file (from get_mirca_read_counts.py)",metavar="file"),
  make_option(c("-o","--outf"),type='character',default="mirca_results.csv.gz",help="output file [%default]",metavar="file"),
  make_option(c("-c","--conditions"),type="character",help="comma-separated list of conditions (corresponding to the columns in input file)",metavar="list"),
  make_option("--reference",type='character',help="reference condition (default: alphabetically first)",metavar="string"),
  make_option("--num_covariate",type="character",help="comma-separated list of numerical covariates",metavar="list"),
  make_option("--cat_covariate",type="character",help="comma-separated list of categorical covariates",metavar="list"),
  make_option("--alpha",type="double",default=0.1,help="significance cutoff alpha [%default]",metavar="double"),
  make_option("--test",type="character",default='LRT',help="statistical test: Wald (with betaPrior) or LRT [%default]",metavar="string")
)

opt_parser <- OptionParser(option_list = option_list)
options <- parse_args(opt_parser)

if (is.null(options$input)) {
  print_help(opt_parser)
  stop("please supply input file",call.=FALSE)
}

if (is.null(options$conditions)) {
  print_help(opt_parser)
  stop("please supply comma-separated list of conditions",call.=FALSE)
}

conditions <- strsplit(options$conditions,',',fixed=TRUE)[[1]]

print(paste("reading from",options$input))
counts.all <- read.csv(options$input,sep='\t',comment='#',header=T)

sampleName <- colnames(counts.all)[4:length(colnames(counts.all))]
if (length(conditions) != length(sampleName)) {
  stop("number of conditions (",length(conditions),") doesn't match number of columns (",length(sampleName),") in ",options$input,call.=FALSE)
}

if (options$test=='Wald') {
  full_design <- ~group
} else {
  full_design <- ~condition + motif + condition:motif
  reduced_design <- ~condition + motif
}

# ADD COVARIATES FOR MORE COMPLEX DESIGN IF AVAILABLE
has_covariate <- FALSE

if (!is.null(options$num_covariate)) {
  has_covariate <- TRUE
  print('using numerical covariate')
  covariate <- as.numeric(strsplit(options$num_covariate,',',fixed=TRUE)[[1]])
} else if (!is.null(options$cat_covariate)) {
  has_covariate <- TRUE
  print('using categorical covariate')
  covariate <- strsplit(options$cat_covariate,',',fixed=TRUE)[[1]]
}

if (has_covariate) {
  if (length(covariate) != length(sampleName)) {
    stop("number of covariates (",length(covariate),") doesn't match number of columns (",length(sampleName),") in ",options$input,call.=FALSE)
  } 
  if (options$test=="Wald") {
    full_design <- ~group + covariate
  } else{
    full_design <- ~condition + motif + covariate + condition:motif
    reduced_design <- ~ condition + covariate + motif
  }
}

# take only named conditions, drop empty ones or those called "_"
take <- (conditions!='_') & (conditions!='')
colData <- data.frame(condition=conditions[take],row.names=sampleName[take])
if (has_covariate) {
   colData['covariate'] <- covariate[take]
}
if (is.null(options$reference)) {
  reference.condition <- sort(unique(conditions[take]))[1]
} else {
  reference.condition <- options$reference
}

columns <- c('gene',sampleName[take])
if (has_covariate) {
   print(paste('using columns',paste(columns[2:length(columns)],collapse=','),'for conditions',paste(conditions[take],collapse=','),'and covariates',paste(covariate[take],collapse=',')))
} else {
   print(paste('using columns',paste(columns[2:length(columns)],collapse=','),'for conditions',paste(conditions[take],collapse=',')))
}

all.motifs <- setdiff(unique(counts.all$motif),'tot')
all.conds <- setdiff(unique(conditions[take]),reference.condition)

counts.tot <- data.table(na.omit(counts.all[counts.all$motif=='tot',columns]))
counts.tot.df <- as.data.frame.matrix(counts.tot[,2:length(counts.tot),with=FALSE],row.names=as.vector(counts.tot$gene))
dds.tot <- DESeqDataSetFromMatrix(counts.tot.df,colData,~condition)
sf.tot <- estimateSizeFactorsForMatrix(counts(dds.tot))

res <- list()
for (mot in all.motifs) {
  print (paste0('deseq2 (',options$test,') for ',mot))
  counts.motif <- data.table(na.omit(counts.all[counts.all$motif==mot,columns]))
  counts.here <- merge(counts.tot,counts.motif,by='gene',all.y=TRUE,suffixes=c('_tot',paste0('_',mot)))
  counts.here.df <- as.data.frame.matrix(counts.here[,2:length(counts.here),with=FALSE],row.names=as.vector(counts.here$gene))
  sampleName <- colnames(counts.here)[2:length(colnames(counts.here))]
  motifs <- c(rep("tot",sum(take)),rep(mot,sum(take)))
  # EDIT COLDATA AND DESIGN IF MORE COVARIATES ARE AVAILABLE
  colData <- data.frame(motif=motifs,condition=c(conditions[take],conditions[take]),
                        group=paste0(c(conditions[take],conditions[take]),motifs),
                        row.names=sampleName)
  if (has_covariate) {
    colData['covariate'] <- c(covariate[take],covariate[take])
  }
  dds <- DESeqDataSetFromMatrix(counts.here.df, colData, full_design)
  dds$motif <- factor(dds$motif, levels=c('tot',mot))
  dds$condition <- factor(dds$condition, levels=c(reference.condition,all.conds))
  # filter out genes with zero counts at this motif
  dds <- dds[rowSums(counts(dds)[,dds$motif==mot]) > 0,]
  sizeFactors(dds) <- c(sf.tot,sf.tot)
  tryCatch({
    if (options$test=="Wald") {
      dds <- DESeq(dds, test='Wald')
    } else {
      dds <- DESeq(dds, test='LRT', reduced=reduced_design)
    }
    for (cond in all.conds) {
      print (paste0('results for ',cond,' vs ',reference.condition))
      if (options$test=="Wald") {
        # use contrast: compare count ratio motif/tot at cond to ratio motif/tot at reference.condition
        # this works for DESeq2 v1.14.1
        # contrast <- list(c(paste0('group',cond,gsub('+','.',mot,fixed=TRUE)),
        #                    paste0('group',reference.condition,'tot')),
        #                  c(paste0('group',reference.condition,gsub('+','.',mot,fixed=TRUE)),
        #                    paste0('group',cond,'tot')))
        # this works for DESeq2 v1.16.1
        contrast <- list(c(paste0('group_',reference.condition,'tot_vs_',reference.condition,gsub('+','.',mot,fixed=TRUE)),
                           paste0('group_',cond,gsub('+','.',mot,fixed=TRUE),'_vs_',reference.condition,gsub('+','.',mot,fixed=TRUE))),
                         c(paste0('group_',cond,'tot_vs_',reference.condition,gsub('+','.',mot,fixed=TRUE))))
        tmp <- as.data.frame(results(dds,contrast=contrast,alpha=options$alpha))
      } else {
        # use interaction term for motif and cond
        tmp <- as.data.frame(results(dds,name=paste0('condition',cond,'.motif',gsub('+','.',mot,fixed=TRUE)),alpha=options$alpha))
      }
      tmp <- data.table(gene=row.names(tmp),tmp)
      setkey(tmp,'gene')
      colnames(tmp) <- c('gene',paste(mot,cond,colnames(tmp)[2:length(tmp)],sep='.'))
      res[[paste0(mot,'.',cond)]] <- tmp
    }
  }, error=function(err) {
    print(paste0('skipping ',mot,': ',err))
  })
}  

print('combining results')
res <- as.data.frame(Reduce(function(...) merge(...,all=T),res))
print('fixing row names')
row.names(res) <- res$gene
print('selecting columns')
res <- res[,2:length(res)]
print(paste('writing combined results to',options$outf))
if (grep("*.gz$",options$outf)) {
  write.csv(res,file=gzfile(options$outf))
} else {
  write.csv(res,file=options$outf)
}

