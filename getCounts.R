#makeComparison.R
library(Rsubread)
library(BiocParallel)


#------------------------------------------------------
n.procs=12
register(MulticoreParam(workers=n.procs))
#------------------     READ ARGS  --------------------
options(echo=TRUE)
#------------------------------------------------------
args=commandArgs(trailingOnly=TRUE)
#------------------------------------------------------
print(args)

BAML=args[1]
ANNO=args[2]

BAM=read.table(BAML)

CNTS=featureCounts(BAM,
                   annot.ext        =ANNO,
                   GTF.featureType  ="unknown",
                   GTF.attrType     ="gene_id",
                   allowMultiOverlap=F,
                   minOverlap       =1,
                   isPairedEnd      =TRUE,
                   countMultiMappingReads=F,
                   strandSpecific   =0,
                   juncCounts       =T,
                   genome           ="mm10",
                   nthreads         =12,
                   PE_orientation   ="fr",
                   reportReads      =T)

save.image(file="~/Projects/counts.Rda")
