library(HTSanalyzeR)
library(KEGG.db)
library(GO.db)
library(BiocParallel)
library(MeSH.db)
#------------------------------------------------------
n.procs=12
register(MulticoreParam(workers=n.procs))

#---------------------------------------------------------------
expr=args[1]
gmtf=args[2]
thresh=args[3]
gsea=args[4]
gsoa=args[5]
org=args[6]
#---------------------------------------------------------------
KEGG=KeggGeneSets(species=org)
GO$mf=GOGeneSets(species=org,"MF")
GO$bp=GOGeneSets(species=org,"BP")

readGMT<-function(gmtfile)
{
  gmt=NULL
  gmt=readLines(gmtfile)
  gmt=lapply(gmt,function(x){return(unlist(strsplit(x,"\t")))})
  nms=lapply(gmt,function(x){return(x[1])})
  gmt=lapply(gmt,function(x){return(tolower(x[3:length(x)]))})
  
  names(gmt)=nms
  return(gmt)
}


gmt=readGMT(gmtf)
dat=read.table(expr,header=T,row.names=1,sep="\t")
logfc=dat$logfc
names(logfc)=tolower(rownames(dat))
#--------------------------------------------------------------------
#
#  reMAP data to KEGG, GOmf GObp to gene symbols
#
#--------------------------------------------------------------------
signif=rownames(expr)[which(expr$pval<thresh)]
allGMT=list(kegg=KEGG,gomf=GO$mf,gobp=GO$bp,msigdb=gmt)
namesGMT=names(allGMT)
results=analyzeGeneSetCollections(geneList=logfc, 
                                  hits=signif, 
                                  minGeneSetSize=5, 
                                  allGMT,
                                  doGSEA=gsea,
                                  doGSOA=gsoa)

gsca=new("GSCA",listOfGeneSetCollections=allGMT,
         para=list(pValueCutoff=0.05,pAdjustMethod="BH",nPermutations=1000,
                   minGeneSetSize=5,exponent=1),
         hits=signif,
         results=results,
         preprocessed=T,
         geneList=logfc)

#---------------------------------------------------------------------
for(i in 1:length(namesGMT))
{
   topHITS=getTopGeneSets(gsca,"GSEA.results",namesGMT[i],allSig=TRUE)
   
   plotGSEA(gsca,gscs=namesGMT[i],
            allSig=TRUE,
            filepath=out,
            output="pdf")
   
}
#----------------------------------------------------------------------

