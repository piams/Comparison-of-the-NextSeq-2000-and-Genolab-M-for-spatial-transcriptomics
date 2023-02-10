## loading libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(writexl))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(venn))
suppressPackageStartupMessages(library(pheatmap))

working_dir<-'/home/piam/analysis/Comparison_10x'
setwd(working_dir)

file_list<-read.table(file =file.path(working_dir,'path.txt'),sep='\t',header=T)

file_list$Seq[grepl('new',file_list$Seq,ignore.case = T)]<-gsub('New','Genolab',file_list$Seq[grepl('new',file_list$Seq,ignore.case = T)])
file_list$Slice[!grepl('new',file_list$Slice,ignore.case = T)]<-paste0(file_list$Slice[!grepl('new',file_list$Slice,ignore.case = T)],'_Illumina')
file_list$Slice[grepl('new',file_list$Slice,ignore.case = T)]<-gsub('new','Genolab',file_list$Slice[grepl('new',file_list$Slice,ignore.case = T)])

#Open 10x visium files

n=1
for (i in file_list$Path){
  assign(file_list$Slice[n],
    Load10X_Spatial(data.dir = i,slice = file_list$Slice[n]))
  n<-n+1
};rm(i);rm(n)

#Add metadata

n=1
for (i in file_list$Slice){
  a<-get(i)
  a@meta.data$group<-file_list$Seq[n]
  a$orig.ident<-i
  a@active.ident<-as.factor(a$orig.ident)
  assign(i,a)
  n<-n+1
};rm(i);rm(n);rm(a)

st_list<-list(A11_Illumina,A12_Illumina,B12_Illumina,A11_Genolab,A12_Genolab,B12_Genolab)
#
cbind_df<-function(x,y){
  colnames(x@meta.data)<-paste0(colnames(x@meta.data),gsub('[A-Z][0-9]{2}','',paste(quote(x))))
  colnames(y@meta.data)<-paste0(colnames(y@meta.data),gsub('[A-Z][0-9]{2}','',paste(quote(y))))
  x<-x@meta.data
  y<-y@meta.data
  x<-cbind(x,y)
  return(x)}

A11<-cbind_df(A11_Illumina,A11_Genolab)
A12<-cbind_df(A12_Illumina,A12_Genolab)
B12<-cbind_df(B12_Illumina,B12_Genolab)

cor(A11$nCount_Spatialx,A11$nCount_Spatialy,method = 'pearson')^2
cor(A12$nCount_Spatialx,A12$nCount_Spatialy,method = 'pearson')^2
cor(B12$nCount_Spatialx,B12$nCount_Spatialy,method = 'pearson')^2

cor(A11$nFeature_Spatialx,A11$nFeature_Spatialy,method = 'pearson')^2
cor(A12$nFeature_Spatialx,A12$nFeature_Spatialy,method = 'pearson')^2
cor(B12$nFeature_Spatialx,B12$nFeature_Spatialy,method = 'pearson')^2

svg('cor_umi_genes.svg',7,5)
ggplot(A11,aes(x=nCount_Spatialx,y=nCount_Spatialy))+geom_point(col=brewer.pal(8,'Set1')[1])+
  xlab('Illumina UMI')+ylab('Genolab UMI')+theme_bw()+labs(title = 'A1-1',subtitle = as.expression(bquote("" ~ R^2 ~ "= 0.987")))+
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
ggplot(A12,aes(x=nCount_Spatialx,y=nCount_Spatialy))+geom_point(col=brewer.pal(8,'Set1')[1])+
  xlab('Illumina UMI')+ylab('Genolab UMI')+theme_bw()+labs(title = 'A1-2',subtitle = as.expression(bquote("" ~ R^2 ~ "= 0.998")))+
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  #theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
ggplot(B12,aes(x=nCount_Spatialx,y=nCount_Spatialy))+geom_point(col=brewer.pal(8,'Set1')[1])+
  xlab('Illumina UMI')+ylab('Genolab UMI')+theme_bw()+labs(title = 'B1-2',subtitle = as.expression(bquote("" ~ R^2 ~ "= 0.987")))+
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  #theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
ggplot(A11,aes(x=nFeature_Spatialx,y=nFeature_Spatialy))+geom_point(col=brewer.pal(8,'Set1')[2])+
  xlab('Illumina Genes')+ylab('Genolab Genes')+theme_bw()+labs(subtitle = as.expression(bquote("" ~ R^2 ~ "= 0.990")))+
  theme(plot.subtitle = element_text(hjust = 0.5))+
ggplot(A12,aes(x=nFeature_Spatialx,y=nFeature_Spatialy))+geom_point(col=brewer.pal(8,'Set1')[2])+
  xlab('Illumina Genes')+ylab('Genolab Genes')+theme_bw()+labs(subtitle = as.expression(bquote("" ~ R^2 ~ "= 0.997")))+
  #theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(plot.subtitle = element_text(hjust = 0.5))+
ggplot(B12,aes(x=nFeature_Spatialx,y=nFeature_Spatialy))+geom_point(col=brewer.pal(8,'Set1')[2])+
  xlab('Illumina Genes')+ylab('Genolab Genes')+theme_bw()+labs(subtitle = as.expression(bquote("" ~ R^2 ~ "= 0.984")))+
  #theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

rbind_df<-function(x,y){
  x<-x@meta.data
  y<-y@meta.data
  x<-rbind(x,y)
  return(x)}

A11<-rbind_df(A11_Illumina,A11_Genolab)
A12<-rbind_df(A12_Illumina,A12_Genolab)
B12<-rbind_df(B12_Illumina,B12_Genolab)

svg('umi_genes_plot.svg',8,3)
ggplot(A11,aes(x=nCount_Spatial,y=nFeature_Spatial,color=group))+geom_point()+scale_color_manual(values=brewer.pal(6,'Set1')[c(2,1)])+
  theme_bw()+xlab('UMI')+ylab('Genes')+labs(title = 'A1-1')+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
ggplot(A12,aes(x=nCount_Spatial,y=nFeature_Spatial,color=group))+geom_point()+scale_color_manual(values=brewer.pal(6,'Set1')[c(2,1)])+
  theme_bw()+xlab('UMI')+ylab('Genes')+labs(title = 'A1-2')+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
ggplot(B12,aes(x=nCount_Spatial,y=nFeature_Spatial,color=group))+geom_point()+scale_color_manual(values=brewer.pal(6,'Set1')[c(2,1)])+
  theme_bw()+xlab('UMI')+ylab('Genes')+labs(title = 'B1-2')+theme(plot.title = element_text(hjust = 0.5))
dev.off()

##Filter genes expressing in less than 10 spots and spots having less than 200 unique genes
st_list_filtered<-mclapply(st_list,function(y){
  counts <- GetAssayData(y, assay = "Spatial")
  nrow(counts)
  counts <- counts[apply(counts,1,function(x) length(which(x>0)))>=10,]
  nrow(counts)
  y<-subset(y,features=rownames(counts))
  length(which((y[['nFeature_Spatial']]>200)==F))
  y<-subset(y,nFeature_Spatial>200)},
  mc.cores = 8)

st_list_filtered<-mclapply(st_list_filtered,function(x) SCTransform(x,variable.features.n = NULL,
                  variable.features.rv.th = 1.1,return.only.var.genes = F,
                  assay = 'Spatial'),mc.cores = 6)
#saveRDS(st_list_filtered,file.path(working_dir,'st_list_filtered.rds'))
#st_list_filtered<-readRDS('st_list_filtered.rds')

lapply(st_list_filtered,function(x){
  name<-x@images%>%names
  df<-data.frame(col=colnames(x))
  colnames(df)<-name
  write.table(df,paste0(name,'.tsv'),sep='\t',row.names = F)
})

  
data.combined<-Reduce(function(x,y) merge(x,y),st_list_filtered)

data.combined<-SCTransform(data.combined,variable.features.n = NULL,
                           variable.features.rv.th = 1.1,return.only.var.genes = F,
                           assay = 'Spatial')

DefaultAssay(data.combined)<-'SCT'

data.combined<-RunPCA(data.combined,npcs = 50)
data.combined<-FindNeighbors(data.combined,dims = 1:30)
data.combined<-FindClusters(data.combined,resolution = 0.5)
data.combined <- RunUMAP(data.combined, dims = 1:30)
data.combined <- RunTSNE(data.combined, dims = 1:30)
#saveRDS(data.combined,file.path(working_dir,'data.combined.rds'))
#data.combined<-readRDS(file.path(working_dir,'data.combined.rds'))

data.combined@meta.data

variable_genes<-data.combined@assays$SCT@var.features

bulk<-bind_cols(lapply(levels(data.combined@meta.data$seurat_clusters),function(x){
I<-data.combined[,rownames(subset(data.combined@meta.data,seurat_clusters==x & group=='Illumina'))]
G<-data.combined[,rownames(subset(data.combined@meta.data,seurat_clusters==x & group=='Genolab'))]
I<-apply(I@assays$SCT@counts,1,sum)
G<-apply(G@assays$SCT@counts,1,sum)
df<-data.frame(I=I,G=G)
colnames(df)<-paste0(x,c('_I','_G'))
return(df)
}))

bulk_var<-bulk[variable_genes,]
bulk_var<-as.data.frame(t(bulk_var))

dist_bulk<-dist(bulk_var)
plot(hclust(dist_bulk))

#pheatmap(cor(bulk_var))
#pheatmap(log10(bulk_var+1))

spots<-lapply(unique(data.combined@meta.data$seurat_clusters),function(x){
  a <- colnames(subset(data.combined,seurat_clusters==x&group=='Illumina'))
  b <- colnames(subset(data.combined,seurat_clusters==x&group=='Genolab'))
  a<-gsub('_.*','',a)
  b<-gsub('_.*','',b)
  list<-list(a,b)
  names(list)<-paste0(x,c('_I','_G'))
  return(list)
})


splited.data <- SplitObject(data.combined, split.by = "orig.ident")


svg(file.path(working_dir,'umap.svg'),10,8)
(DimPlot(subset(data.combined,group=='Illumina'), reduction = "umap",label = T,raster.dpi = c(2024,2024),pt.size = 9,
         raster = T,cols = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2')))+
  DimPlot(data.combined, reduction = "umap", group.by = c("orig.ident"),label = F,raster.dpi = c(2024,2024),pt.size = 9,
          raster = T,cols = brewer.pal(6,'Set1'))) /
(DimPlot(subset(data.combined,group=='Genolab'), reduction = "umap",label = T,raster.dpi = c(2024,2024),pt.size = 9,
         raster = T,cols = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2')))+
  DimPlot(data.combined, reduction = "umap", group.by = c("group"),label = F,raster.dpi = c(2024,2024),pt.size = 9,
          raster = T,cols = brewer.pal(6,'Set1')[2:1]))
dev.off()





plot3<-DimPlot(splited.data[[3]], reduction = "umap", group.by = c("ident"),label = T,
               raster = T,pt.size = 15,raster.dpi = c(2024,2024),
               cols = c(brewer.pal(9,'Set1')[c(2,6,9)],brewer.pal(8,'Dark2')[c(4,7)]))/
  DimPlot(splited.data[[6]], reduction = "umap", group.by = c("ident"),label = T,
          raster = T,pt.size = 15,raster.dpi = c(2024,2024),
          cols = c(brewer.pal(9,'Set1')[c(2,6,9)],brewer.pal(8,'Dark2')[c(4,7)]))

plot2<-DimPlot(splited.data[[2]], reduction = "umap", group.by = c("ident"),label = T,
               raster = T,pt.size = 15,raster.dpi = c(2024,2024),
               cols = c(brewer.pal(9,'Set1')[c(1,5)],brewer.pal(8,'Dark2')[c(1,3)]))/
  DimPlot(splited.data[[5]], reduction = "umap", group.by = c("ident"),label = T,
          raster = T,pt.size = 15,raster.dpi = c(2024,2024),
          cols = c(brewer.pal(9,'Set1')[c(1,5)],brewer.pal(8,'Dark2')[c(1,3)]))

plot1<-DimPlot(splited.data[[1]], reduction = "umap", group.by = c("ident"),label = T,
               raster = T,pt.size = 15,raster.dpi = c(2024,2024),
               cols = c(brewer.pal(9,'Set1')[c(3,4,7,8)],brewer.pal(8,'Dark2')[c(2,5,6)]))/
  DimPlot(splited.data[[4]], reduction = "umap", group.by = c("ident"),label = T,
          raster = T,pt.size = 15,raster.dpi = c(2024,2024),
          cols = c(brewer.pal(9,'Set1')[c(3,4,7,8)],brewer.pal(8,'Dark2')[c(2,5,6)]))

svg(file.path(working_dir,'single_umap.svg'),10,5)
plot1|plot2|plot3
dev.off()





svg(filename = file.path(working_dir,'venn_genes.svg'),width = 10,height = 3)
venn(x=list(rownames(st_list_filtered[[1]]),rownames(st_list_filtered[[4]])),ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-1",size=4)+
venn(x=list(rownames(st_list_filtered[[2]]),rownames(st_list_filtered[[5]])),ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-2",size=4)+
venn(x=list(rownames(st_list_filtered[[3]]),rownames(st_list_filtered[[6]])),ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "B1-2",size=4)
dev.off()


A11_Genolab<-st_list_filtered[[4]][which(rownames(st_list_filtered[[4]]) %in% rownames(st_list_filtered[[1]])),
                                   which(colnames(st_list_filtered[[4]]) %in% colnames(st_list_filtered[[1]]))]
A12_Genolab<-st_list_filtered[[5]][which(rownames(st_list_filtered[[5]]) %in% rownames(st_list_filtered[[2]]))
                                   ,which(colnames(st_list_filtered[[5]]) %in% colnames(st_list_filtered[[2]]))]
B12_Genolab<-st_list_filtered[[6]][which(rownames(st_list_filtered[[6]]) %in% rownames(st_list_filtered[[3]]))
                                   ,which(colnames(st_list_filtered[[6]]) %in% colnames(st_list_filtered[[3]]))]
A11_Illumina<-st_list_filtered[[1]][which(rownames(st_list_filtered[[1]]) %in% rownames(st_list_filtered[[4]])),]
A12_Illumina<-st_list_filtered[[2]][which(rownames(st_list_filtered[[2]]) %in% rownames(st_list_filtered[[5]])),]
B12_Illumina<-st_list_filtered[[3]][which(rownames(st_list_filtered[[3]]) %in% rownames(st_list_filtered[[6]])),]

n=1;rA11<-c()
while (n<=ncol(A11_Genolab)) {
  rA11<<-c(rA11,cor(A11_Genolab@assays$SCT@counts[,n],A11_Illumina@assays$SCT@counts[,n]))
  n<-n+1
}
n=1;rA12<-c()
while (n<=ncol(A12_Genolab)) {
  rA12<<-c(rA12,cor(A12_Genolab@assays$SCT@counts[,n],A12_Illumina@assays$SCT@counts[,n]))
  n<-n+1
}
n=1;rB12<-c()
while (n<=ncol(B12_Genolab)) {
  rB12<<-c(rB12,cor(B12_Genolab@assays$SCT@counts[,n],B12_Illumina@assays$SCT@counts[,n]))
  n<-n+1
}

rA11<-data.frame(rA11=rA11)
rA12<-data.frame(rA12=rA12)
rB12<-data.frame(rB12=rB12)

mean(rA11$rA11)
mean(rA12$rA12)
mean(rB12$rB12)


svg('spot_cor.svg',8,3)
ggplot(rA11,aes(x=rA11))+geom_histogram(binwidth = 0.005,fill=brewer.pal(6,'Set1')[1],alpha = 0.8)+
  xlab('Pearson correlation')+ylab('Number of spots')+theme_bw()+labs(title = 'A1-1')+
  theme(plot.title = element_text(hjust = 0.5))+
ggplot(rA12,aes(x=rA12))+geom_histogram(binwidth = 0.004,fill=brewer.pal(6,'Set1')[1],alpha = 0.8)+
  xlab('Pearson correlation')+ylab('Number of spots')+theme_bw()+labs(title = 'A1-2')+
  theme(plot.title = element_text(hjust = 0.5))+
ggplot(rB12,aes(x=rB12))+geom_histogram(binwidth = 0.001,fill=brewer.pal(6,'Set1')[1],alpha = 0.8)+
  xlab('Pearson correlation')+ylab('Number of spots')+theme_bw()+labs(title = 'B1-2')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



DE<-mclapply(splited.data,function(y){
  x<-FindAllMarkers(y,assay = 'SCT',logfc.threshold = 0.25,return.thresh = 0.01,only.pos = T)
  x<-subset(x,p_val_adj<0.01)
  rownames(x)<-NULL
  x$group<-rep(unique(y@meta.data$group),nrow(x))
  return(x)
},mc.cores = 30)
#saveRDS(DE,'DE.rds')
#DE<-readRDS('DE.rds')



col<-c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'))
names(col)<-levels(data.combined$seurat_clusters)

svg('slide_cols.svg',20,8)
SpatialDimPlot(data.combined,cols = col)
dev.off()


slide_degs<-sapply(names(DE)[1:3],function(x){
  x<-gsub('_.*$','',x)
  x<-DE[grep(x,names(DE))]
  x<-rbind(x[[1]],x[[2]])
  n<-levels(x$cluster)
  n<-sapply(n,function(y){
    x<-subset(x,cluster==y)
    x<-list(subset(x,group=='Illumina')$gene,subset(x,group=='Genolab')$gene)
    return(x)
    })
  return(n)
  })

#Venn between clusters in A11
svg('cluster_venn.svg',10,10)
venn(x=slide_degs[["A11_Illumina"]][,c('2')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-1\n cluster 2",size=4)+
venn(x=slide_degs[["A11_Illumina"]][,c('3')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-1\n cluster 3",size=4)+
venn(x=slide_degs[["A11_Illumina"]][,c('6')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-1\n cluster 6",size=4)+
venn(x=slide_degs[["A11_Illumina"]][,c('7')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-1\n cluster 7",size=4)+
venn(x=slide_degs[["A11_Illumina"]][,c('10')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-1\n cluster 10",size=4)+
venn(x=slide_degs[["A11_Illumina"]][,c('13')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-1\n cluster 13",size=4)+
venn(x=slide_degs[["A11_Illumina"]][,c('14')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-1\n cluster 14",size=4)+

#Venn between clusters in A12

venn(x=slide_degs[["A12_Illumina"]][,c('0')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-2\n cluster 0",size=4)+
venn(x=slide_degs[["A12_Illumina"]][,c('4')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-2\n cluster 4",size=4)+
venn(x=slide_degs[["A12_Illumina"]][,c('9')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-2\n cluster 9",size=4)+
venn(x=slide_degs[["A12_Illumina"]][,c('11')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "A1-2\n cluster 11",size=4)+

  
#Venn between clusters in B12

venn(x=slide_degs[["B12_Illumina"]][,c('1')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "B1-2\n cluster 1",size=4)+
venn(x=slide_degs[["B12_Illumina"]][,c('5')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "B1-2\n cluster 5",size=4)+
venn(x=slide_degs[["B12_Illumina"]][,c('8')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "B1-2\n cluster 8",size=4)+
venn(x=slide_degs[["B12_Illumina"]][,c('12')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "B1-2\n cluster 12",size=4)+
venn(x=slide_degs[["B12_Illumina"]][,c('15')],ilcs = 0.9,sncs = 0.9,plotsize = 15,
  snames = c('Illumina','Genolab'),box = F,zcolor = brewer.pal(6,'Set1')[1:2],opacity = 0.6,ggplot = T)+theme_void()+
  annotate("text", x = 490, y = 870, label = "B1-2\n cluster 15",size=4)
dev.off()



genes<-lapply(slide_degs,function(x){
  list<-lapply(colnames(x),function(c){
    y<-x[,c]
    counts_illumina<-apply(subset(data.combined,seurat_clusters==c&group=="Illumina")@assays$Spatial@counts,1,sum)
    counts_genolab<-apply(subset(data.combined,seurat_clusters==c&group=="Genolab")@assays$Spatial@counts,1,sum)
    
    de<-rbind(DE$A11_Illumina,DE$A12_Illumina,DE$B12_Illumina,DE$A11_Genolab,DE$A12_Genolab,DE$B12_Genolab)
    de_illumina<-subset(de,cluster==c&group=="Illumina")
    de_genolab<-subset(de,cluster==c&group=="Genolab")
    
    #Common illumina counts
    common_illumina<-data.frame(counts=counts_illumina[y[[1]][duplicated(c(y[[1]],y[[2]]),fromLast = T)]])
    common_illumina$group<-rep('Common Illumina',nrow(common_illumina))
    common_illumina$gene<-y[[1]][duplicated(c(y[[1]],y[[2]]),fromLast = T)]
    common_illumina$FDR<-de_illumina[which(de_illumina$gene %in% rownames(common_illumina)),]$p_val_adj
    common_illumina$LFC<-de_illumina[which(de_illumina$gene %in% rownames(common_illumina)),]$avg_log2FC
    rownames(common_illumina)<-NULL
    
    #Common genolab counts
    common_genolab<-data.frame(counts=counts_genolab[y[[2]][duplicated(c(y[[2]],y[[1]]),fromLast = T)]])
    common_genolab$group<-rep('Common Genolab',nrow(common_genolab))
    common_genolab$gene<-y[[2]][duplicated(c(y[[2]],y[[1]]),fromLast = T)]
    common_genolab$FDR<-de_genolab[which(de_genolab$gene %in% rownames(common_genolab)),]$p_val_adj
    common_genolab$LFC<-de_genolab[which(de_genolab$gene %in% rownames(common_genolab)),]$avg_log2FC
    rownames(common_genolab)<-NULL
    
    illumina<-data.frame(counts=counts_illumina[y[[1]][-which(y[[1]] %in% y[[1]][duplicated(c(y[[1]],y[[2]]),fromLast = T)])]])
    illumina$group<-rep('Unique Illumina',nrow(illumina))
    illumina$gene<-y[[1]][-which(y[[1]] %in% y[[1]][duplicated(c(y[[1]],y[[2]]),fromLast = T)])]
    illumina$FDR<-de_illumina[which(de_illumina$gene %in% rownames(illumina)),]$p_val_adj
    illumina$LFC<-de_illumina[which(de_illumina$gene %in% rownames(illumina)),]$avg_log2FC
    rownames(illumina)<-NULL
    
    genolab<-data.frame(counts=counts_genolab[y[[2]][-which(y[[2]] %in% y[[1]][duplicated(c(y[[1]],y[[2]]),fromLast = T)])]])
    genolab$group<-rep('Unique Genolab',nrow(genolab))
    genolab$gene<-y[[2]][-which(y[[2]] %in% y[[1]][duplicated(c(y[[1]],y[[2]]),fromLast = T)])]
    genolab$FDR<-de_genolab[which(de_genolab$gene %in% rownames(genolab)),]$p_val_adj
    genolab$LFC<-de_genolab[which(de_genolab$gene %in% rownames(genolab)),]$avg_log2FC
    rownames(genolab)<-NULL
    
    y<-rbind(common_illumina,common_genolab,illumina,genolab)
    return(y)})
  names(list)<-colnames(x)
  return(list)
})


#A11 overlapped and unique genes counts
ggplot(genes$A11_Illumina$'2',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 2')+
ggplot(genes$A11_Illumina$'3',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 3')+
ggplot(genes$A11_Illumina$'6',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 6')+
ggplot(genes$A11_Illumina$'7',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 7')+
ggplot(genes$A11_Illumina$'10',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 10')+
ggplot(genes$A11_Illumina$'13',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 13')+
ggplot(genes$A11_Illumina$'14',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 14')

#A12 overlapped and unique genes counts
ggplot(genes$A12_Illumina$'0',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-2 cluster 0')+
ggplot(genes$A12_Illumina$'4',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-2 cluster 4')+
ggplot(genes$A12_Illumina$'9',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-2 cluster 9')+
ggplot(genes$A12_Illumina$'11',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-2 cluster 11')

#B12 overlapped and unique genes counts
ggplot(genes$B12_Illumina$'1',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'B1-2 cluster 1')+
ggplot(genes$B12_Illumina$'5',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'B1-2 cluster 5')+
ggplot(genes$B12_Illumina$'8',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'B1-2 cluster 8')+
ggplot(genes$B12_Illumina$'12',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'B1-2 cluster 12')+
ggplot(genes$B12_Illumina$'15',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.2,alpha=0.8)+theme_bw()+
scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'B1-2 cluster 15')


#A11
svg('cluster_gene_plots.svg',30,15)
ggplot(genes$A11_Illumina$'2',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.05,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
ggplot(genes$A11_Illumina$'2',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=3.75,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 2')+
ggplot(genes$A11_Illumina$'2',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.075,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$A11_Illumina$'3',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.05,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$A11_Illumina$'3',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=10,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 3')+
  ggplot(genes$A11_Illumina$'3',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.085,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$A11_Illumina$'6',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.075,alpha=0.9)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$A11_Illumina$'6',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=10,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 6')+
  ggplot(genes$A11_Illumina$'6',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.125,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$A11_Illumina$'7',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.03,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$A11_Illumina$'7',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=4,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 7')+
  ggplot(genes$A11_Illumina$'7',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.1,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$A11_Illumina$'10',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.04,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$A11_Illumina$'10',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=5,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 10')+
  ggplot(genes$A11_Illumina$'10',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.09,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$A11_Illumina$'13',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.125,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$A11_Illumina$'13',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=4,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 13')+
  ggplot(genes$A11_Illumina$'13',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.1,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$A11_Illumina$'14',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.075,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$A11_Illumina$'14',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=2.5,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-1 cluster 14')+
  ggplot(genes$A11_Illumina$'14',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.085,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+


#A12

ggplot(genes$A12_Illumina$'0',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.0145,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$A12_Illumina$'0',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=5.4,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-2 cluster 0')+
  ggplot(genes$A12_Illumina$'0',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.175,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$A12_Illumina$'4',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.038,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$A12_Illumina$'4',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=5.395,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-2 cluster 4')+
  ggplot(genes$A12_Illumina$'4',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.1,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$A12_Illumina$'9',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.05,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$A12_Illumina$'9',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=5,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-2 cluster 9')+
  ggplot(genes$A12_Illumina$'9',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.138,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$A12_Illumina$'11',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.085,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$A12_Illumina$'11',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=5.13,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'A1-2 cluster 11')+
  ggplot(genes$A12_Illumina$'11',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.1,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  

#B12

ggplot(genes$B12_Illumina$'1',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.0112,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$B12_Illumina$'1',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=3.2,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'B1-2 cluster 1')+
  ggplot(genes$B12_Illumina$'1',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.097,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$B12_Illumina$'5',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.065,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$B12_Illumina$'5',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=10,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'B1-2 cluster 5')+
  ggplot(genes$B12_Illumina$'5',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.13,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$B12_Illumina$'8',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.045,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$B12_Illumina$'8',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=5,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'B1-2 cluster 8')+
  ggplot(genes$B12_Illumina$'8',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.075,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$B12_Illumina$'12',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.085,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$B12_Illumina$'12',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=4.48,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'B1-2 cluster 12')+
  ggplot(genes$B12_Illumina$'12',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.153,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  
  ggplot(genes$B12_Illumina$'15',aes(x=LFC,fill=group))+geom_histogram(binwidth=0.11,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+
  ggplot(genes$B12_Illumina$'15',aes(x=-log10(FDR),fill=group))+geom_histogram(binwidth=1.5,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])+labs(title = 'B1-2 cluster 15')+
  ggplot(genes$B12_Illumina$'15',aes(x=log10(counts),fill=group))+geom_histogram(binwidth=0.125,alpha=0.8)+theme_bw()+
  scale_fill_manual(values=brewer.pal(6,'Set1')[1:4])

dev.off()

A11<-genes$A11_Illumina%>%bind_rows()
A12<-genes$A12_Illumina%>%bind_rows()
B12<-genes$B12_Illumina%>%bind_rows()


A11$gene_group<-gsub(' .*','',A11$group)
A12$gene_group<-gsub(' .*','',A12$group)
B12$gene_group<-gsub(' .*','',B12$group)

plot1<-ggplot(A11,aes(x=gene_group,y=-log(FDR),fill=gene_group))+geom_boxplot()+ 
  ggtitle("A11")+xlab('Gene overlapping')+theme_bw()+
  ggplot(A11,aes(x=gene_group,y=LFC,fill=gene_group))+geom_boxplot()+
  ggtitle("A11")+xlab('Gene overlapping')+theme_bw()

plot2<-ggplot(A12,aes(x=gene_group,y=-log(FDR),fill=gene_group))+geom_boxplot()+ 
  ggtitle("A12")+xlab('Gene overlapping')+theme_bw()+
  ggplot(A12,aes(x=gene_group,y=LFC,fill=gene_group))+geom_boxplot()+
  ggtitle("A12")+xlab('Gene overlapping')+theme_bw()

plot3<-ggplot(B12,aes(x=gene_group,y=-log(FDR),fill=gene_group))+geom_boxplot()+ 
  ggtitle("B12")+xlab('Gene overlapping')+theme_bw()+
  ggplot(B12,aes(x=gene_group,y=LFC,fill=gene_group))+geom_boxplot()+
  ggtitle("B12")+xlab('Gene overlapping')+theme_bw()

svg('fdr_lfc.svg',10,8)
plot1/plot2/plot3
dev.off()

#Tables with LFC and FDR distribution btw Common and Unique genes
#LFC
LFC<-list(A11 %>%group_by(gene_group) %>%
summarize(q25 = quantile(LFC, probs = .25), 
            q50 = quantile(LFC, probs = .5),
            q75 = quantile(LFC, probs = .75),
            mean=mean(LFC))%>%as.data.frame(),
A12 %>%group_by(gene_group) %>%
  summarize(q25 = quantile(LFC, probs = .25), 
            q50 = quantile(LFC, probs = .5),
            q75 = quantile(LFC, probs = .75),
            mean=mean(LFC))%>%as.data.frame(),
B12 %>%group_by(gene_group) %>%
  summarize(q25 = quantile(LFC, probs = .25), 
            q50 = quantile(LFC, probs = .5),
            q75 = quantile(LFC, probs = .75),
            mean=mean(LFC))%>%as.data.frame())
#FDR
FDR<-list(A11 %>%group_by(gene_group) %>%
  summarize(q25 = quantile(-log10(FDR), probs = .25), 
            q50 = quantile(-log10(FDR), probs = .5),
            q75 = quantile(-log10(FDR), probs = .75),
            mean=mean(-log10(FDR)))%>%as.data.frame(),
A12 %>%group_by(gene_group) %>%
  summarize(q25 = quantile(-log10(FDR), probs = .25), 
            q50 = quantile(-log10(FDR), probs = .5),
            q75 = quantile(-log10(FDR), probs = .75),
            mean=mean(-log10(FDR)))%>%as.data.frame(),
B12 %>%group_by(gene_group) %>%
  summarize(q25 = quantile(-log10(FDR), probs = .25), 
            q50 = quantile(-log10(FDR), probs = .5),
            q75 = quantile(-log10(FDR), probs = .75),
            mean=mean(-log10(FDR)))%>%as.data.frame())


names(LFC)<-c('A11','A12','B12')
names(FDR)<-c('A11','A12','B12')

write_xlsx(LFC,'LFC.xlsx')
write_xlsx(FDR,'FDR.xlsx')




