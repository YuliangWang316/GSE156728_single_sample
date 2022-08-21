library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
library(tidyverse)
a<-list.files("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
b<-readRDS(a[34])
count_B<-b@assays$data$norm_exprs
count_B<-as.matrix(count_B)
count_B<-as.data.frame(count_B)
count_B_metadata<-as.data.frame(b@colData@listData)
count_B<-count_B[,count_B_metadata$cellID]
colnames(count_B)<-count_B_metadata$cellID.uniq
genename<-rownames(count_B)
d<-str_sub(genename[1],1,4)
if(d == "ENSG"){
  e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
  e<-e[!duplicated(e$SYMBOL),]
  count_B<-count_B[e$ENSEMBL,]
  rownames(count_B)<-e$SYMBOL
}
f<-sort(genename)
if(f[1] == "1"){
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  g<-g[!duplicated(g$SYMBOL),]
  count_B<-count_B[g$ENTREZID,]
  rownames(count_B)<-g$SYMBOL
}
remove(b,e,g,d,f,genename)
e<-readRDS("C:/Users/xjmik/Downloads/data.expression/data/expression/CD4/integration/int.CD4.S35.meta.tb.rds")
rownames(count_B_metadata)<-count_B_metadata$cellID.uniq
h<-intersect(e$cellID.uniq,count_B_metadata$cellID.uniq)
e<-as.data.frame(e)
rownames(e)<-e$cellID.uniq
e_new<-e[h,]
count_B_metadata_new<-count_B_metadata[h,]
metadata<-cbind(e_new,count_B_metadata_new)
remove(e,e_new,count_B_metadata,count_B_metadata_new,h)        #
rownames(metadata)<-colnames(count_B)

# b<-readRDS(a[7])
# count_C<-b@assays$data$norm_exprs
# count_C<-as.matrix(count_C)
# count_C<-as.data.frame(count_C)
# count_C_metadata<-as.data.frame(b@colData@listData)
# count_C<-count_C[,count_C_metadata$cellID]
# colnames(count_C)<-count_C_metadata$cellID.uniq
# genename<-rownames(count_C)
# d<-str_sub(genename[1],1,4)
# if(d == "ENSG"){
#   e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
#   e<-e[!duplicated(e$SYMBOL),]
#   count_C<-count_C[e$ENSEMBL,]
#   rownames(count_C)<-e$SYMBOL
# }
# f<-sort(genename)
# if(f[1] == "1"){
#   g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
#   g<-g[!duplicated(g$SYMBOL),]
#   count_C<-count_C[g$ENTREZID,]
#   rownames(count_C)<-g$SYMBOL
# }
# remove(b,e,g,d,f,genename)
# e<-readRDS("C:/Users/xjmik/Downloads/data.expression/data/expression/CD4/integration/int.CD4.S35.meta.tb.rds")
# rownames(count_C_metadata)<-count_C_metadata$cellID.uniq
# h<-intersect(e$cellID.uniq,count_C_metadata$cellID.uniq)
# e<-as.data.frame(e)
# rownames(e)<-e$cellID.uniq
# e_new<-e[h,]
# count_C_metadata_new<-count_C_metadata[h,]
# metadata_C<-cbind(e_new,count_C_metadata_new)
# remove(e,e_new,count_C_metadata,count_C_metadata_new,h)        #
# rownames(metadata_C)<-colnames(count_C)
# metadata_Total<-rbind(metadata,metadata_C)
# remove(metadata,metadata_C)
# b<-intersect(rownames(count_B),rownames(count_C))
# count_B_new<-count_B[b,]
# count_C_new<-count_C[b,]
# count<-cbind(count_B_new,count_C_new)
# remove(count_B,count_B_new,count_C,count_C_new)
pbmc<-CreateSeuratObject(counts = count_B,meta.data = metadata)
Idents(pbmc)<-pbmc@meta.data$loc
# remove(count,metadata_Total)
# remove(b)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
Idents(pbmc)<-pbmc@meta.data$meta.cluster
remove(all.genes)
b<-levels(Idents(pbmc))
c<-b[which(b == "CD4.c18.Treg.RTKN2" | b == "CD4.c19.Treg.S1PR1" | b == "CD4.c20.Treg.TNFRSF9" | b == "CD4.c21.Treg.ISG")]
Treg<-subset(pbmc,ident=c)
Idents(Treg)<-Treg@meta.data$loc
VlnPlot(Treg,features = "JMJD1C",pt.size = 0,sort = TRUE)


Treg_TP<-subset(Treg,idents = c("T","P"))
VlnPlot(Treg_TP,features = "JMJD1C",pt.size = 0,sort = TRUE)+stat_compare_means()
Treg_TP<-subset(Treg,idents = c("T","N"))
VlnPlot(Treg_TP,features = "JMJD1C",pt.size = 0,sort = TRUE)+stat_compare_means()

Treg<-FindVariableFeatures(Treg, selection.method = "vst", nfeatures = 2000)
Treg<-RunPCA(Treg,features = VariableFeatures(object = Treg))
Treg<-RunUMAP(Treg, dims = 1:20)
Treg<-RunTSNE(Treg, dims = 1:20)

FeaturePlot(Treg,features = "JMJD1C",reduction = "umap",split.by = "loc",order = TRUE) + theme(legend.position = "right")
FeaturePlot(Treg,features = "JMJD1C",reduction = "tsne",split.by = "loc",order = TRUE) + theme(legend.position = "right")

JMJD1C<-FetchData(Treg,vars = "JMJD1C")
JMJD1C_metadata<-Treg@meta.data[rownames(JMJD1C),]
JMJD1C_new<-cbind(JMJD1C,JMJD1C_metadata$loc)
# noise<-rnorm(n=length(x=JMJD1C[,1]))/1e+05
# JMJD1C[,1]<-JMJD1C[,1]+noise
# JMJD1C_new_new<-cbind(JMJD1C,JMJD1C_metadata$loc)
write.table(JMJD1C_new,file = "c:/Users/xjmik/Desktop/TN_34.txt",sep = "\t")
# write.table(JMJD1C_new_new,file = "c:/Users/xjmik/Desktop/BRCATotal.txt",sep = "\t")

Treg_TP<-FindVariableFeatures(Treg_TP, selection.method = "vst", nfeatures = 2000)
Treg_TP<-RunPCA(Treg_TP,features = VariableFeatures(object = Treg_TP))
Treg_TP<-RunUMAP(Treg_TP, dims = 1:20)
Treg_TP<-RunTSNE(Treg_TP, dims = 1:20)

FeaturePlot(Treg_TP,features = "JMJD1C",reduction = "umap",split.by = "loc",order = TRUE) + theme(legend.position = "right")
FeaturePlot(Treg_TP,features = "JMJD1C",reduction = "tsne",split.by = "loc",order = TRUE) + theme(legend.position = "right")

JMJD1C<-FetchData(Treg_TP,vars = "JMJD1C")
JMJD1C_metadata<-Treg_TP@meta.data[rownames(JMJD1C),]
JMJD1C_new<-cbind(JMJD1C,JMJD1C_metadata$loc)
# noise<-rnorm(n=length(x=JMJD1C[,1]))/1e+05
# JMJD1C[,1]<-JMJD1C[,1]+noise
# JMJD1C_new_new<-cbind(JMJD1C,JMJD1C_metadata$loc)
write.table(JMJD1C_new,file = "c:/Users/xjmik/Desktop/TP_20.txt",sep = "\t")
# write.table(JMJD1C_new_new,file = "c:/Users/xjmik/Desktop/BRCATotal.txt",sep = "\t")
