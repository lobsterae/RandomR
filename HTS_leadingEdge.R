#----- leading edge analysis.R ------------------------#
#-------------------------------------------------------
library("DOSE")
library("GOSemSim")
library("ReactomePA")
library("qusage")
#-------------------------------------------------------
data(geneList)
#-------------------------------------------------------
args   =options(trailing=TRUE)
geneSet=read.table(args[1])
org    =args[2]
pval   =args[3]
#-------------------------------------------------------
#--- special case for cancer ---------------------------
geneList=names(geneSet$Gene[which(geneSet$pval < pval)])
universe=geneSet$Gene

gmtfile <- system.file("extdata", "custom.gmt")
c5 <- read.gmt(gmtfile)
c5.parse=NULL

for(i in 1:length(c5)){c5.parse=rbind(c5.parse,cbind(names(c5)[i],c5[[i]]))}

egmt <- enricher(gene, TERM2GENE=c5.parse)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5.parse, verbose=FALSE)
head(egmt2)
dotPlot(egmt2)
#-------------------------------------------------------
#  only for human
#-------------------------------------------------------
CANCER=enrichNCG(names(geneList)     , 
                 pvalueCutoff  = 0.05, 
                 pAdjustMethod = "BH", 
#                universe      =     , 
                 minGSSize     =    5, 
                 qvalueCutoff  =  0.2, 
                 readable      = FALSE)
write.table(CANCER@result,file="CANCER.RESULT")
#-------------------------------------------------------
#  for all organisms
#-------------------------------------------------------
GSEA=gseAnalyzer(geneList, 
                 setType       =    "GO",
                 organism      = "human",
                 nPerm         =    1000, 
                 minGSSize     =       5,
                 pvalueCutoff  =    0.05,
                 pAdjustMethod =    "BH",
                 verbose       =    TRUE)
#-------------------------------------------------------
write.table(GSEA@result,file="GSEA.RESULT")
#-------------------------------------------------------
# for REACTOME
#-------------------------------------------------------
REACT=gsePathway(signif,
                 organism="human",
                 pvalueCutoff=0.05,
                 pAdjustMethod="BH")
write.table(REACT@result,file="REACTOME.RESULT")
#-------------------------------------------------------
# GSEA.SIM.MATRIX=DOSE::doSim(GSEA,GSEA)
#-------------------------------------------------------
# GO.SIM.MATRIX  =GOSemSim::goSim(TERMS1,TERMS2)
#-------------------------------------------------------
# GO.ENRICH      =enrichMap(GSEA,n=50,fixed=TRUE,vertex.label.font=6)
for(i in 1:length(GSEA))
{
  gseaplot(GSEA,   i, by ="all")
}

for(i in 1:length(REACT))
{
  gseaplot(REACT,  i, by ="all")
}

for(i in 1:length(CANCER))
{
 gseaplot(CANCER, i, by ="all")
}
#-------------------------------------------------------
save.image(file=paste(geneSet,".Rda",sep=""))
#-------------------------------------------------------