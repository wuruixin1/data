if (T) {
  dir.create("scripts")
  dir.create("PDFs")
  dir.create("files")
  dir.create("origin_datas/TCGA",recursive = T)
  dir.create("origin_datas/scRNA",recursive = T)
  dir.create("results/files",recursive = T)
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("results/pdf",recursive = T)
  dir.create("results/ana",recursive = T)
}
options(stringsAsFactors = F,check.bounds = F)
#1、单细胞聚类降维和细胞定义#####
dir.create('results/ana/01.scRNA')
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)

dir_name=list.files('origin_datas/scRNA/GSE117570_RAW/',pattern = '.txt.gz')
datalist=list()
for (i in 1:length(dir_name)){
  files = paste0("origin_datas/scRNA/GSE117570_RAW/",dir_name[i])
  counts=fread(file = files,data.table = T,sep = '\t',check.names = F)
  counts=data.frame(counts)
  rownames(counts)=counts[,1]
  counts=counts[,-1]
  rownames(counts) <- gsub("_","-", rownames(counts))
  Samples1=stringr::str_split_fixed(dir_name[i],'_',2)[,1]
  Patient=stringr::str_split_fixed(dir_name[i],'_',3)[,2]
  Type=stringr::str_split_fixed(dir_name[i],'_',4)[,3]
  if(Patient=='P4'){
    cancer='LUSC'
  }else{
    cancer='LUAD'
  }
  colnames(counts)=paste0(Samples1,'_',colnames(counts))
  datalist[[i]]<- CreateSeuratObject(counts=counts,project = Samples1,min.cells = 3, min.features = 250) 
  datalist[[i]] <- AddMetaData(datalist[[i]] , Samples1,col.name = "Samples")
  datalist[[i]] <- AddMetaData(datalist[[i]] , Patient,col.name = "Patient")
  datalist[[i]] <- AddMetaData(datalist[[i]] , Type,col.name = "Type")
  datalist[[i]] <- AddMetaData(datalist[[i]] , cancer,col.name = "cancer")
  
}
names(datalist)=dir_name

for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比
  datalist[[i]] <- sce
  rm(sce)
}

sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])

#细胞数统计
raw_meta=sce@meta.data
raw_count <- table(raw_meta$Samples)
raw_count
sum(raw_count)#8030
pearplot_befor<-VlnPlot(sce,group.by ='Samples', 
                        features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        pt.size = 0, 
                        ncol = 3)
pearplot_befor

#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset =  nFeature_RNA < 5000 & 
              nCount_RNA > 100 &
              percent.mt<15)
})


#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_meta=sce@meta.data

clean_count <- table(clean_meta$Samples)
clean_count
sum(clean_count)#7484
pearplot_after <- VlnPlot(sce,group.by ='Samples', 
                          features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                          pt.size = 0, 
                          ncol = 3)
pearplot_after
save(datalist,file = 'origin_datas/scRNA/datalist.RData')

summary_cells <- as.data.frame(cbind(raw_count,clean_count))
write.table(summary_cells,'results/ana/01.scRNA/cell_count.txt',quote = F,row.names = T,sep='\t')
#降维
load('origin_datas/scRNA/datalist.RData')
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])

sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000,
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
sce <- ScaleData(sce, features =  rownames(sce))
sce <- RunPCA(sce, features = VariableFeatures(sce)) 

dimplot <- DimPlot(sce, reduction = "pca",group.by = 'Samples') 
elbowplot <- ElbowPlot(sce, ndims=50, reduction="pca") 
sc_pca <- dimplot+elbowplot
sc_pca
Dims <- 30
#sce <- RunUMAP(sce, dims=1:Dims, reduction="pca")
sce <- RunTSNE(sce, dims=1:Dims, reduction="pca")

ump.raw<-DimPlot(sce,group.by='Samples',
                 reduction="tsne",
                 label = "T", 
                 pt.size = 0.2,
                 label.size = 0)+
  ggtitle('Before batch elimination')+
  theme(legend.position = 'none')

ump.raw
#ggsave('results/ana/01.scRNA/ump.raw.pdf',ump.raw,height = 7,width = 7)
#去批次
# Normalizing the data
load('origin_datas/scRNA/datalist.RData')
for (i in 1:length(datalist)){
  datalist[[i]]<-NormalizeData(datalist[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  datalist[[i]]<-FindVariableFeatures(datalist[[i]], 
                                      selection.method = "vst", 
                                      nfeatures = 2000,
                                      mean.cutoff=c(0.0125,3),
                                      dispersion.cutoff =c(1.5,Inf))
}

#使用CCA的方法进行剔除批次效应
datalist <- FindIntegrationAnchors(object.list = datalist, dims = 1:40,
                                   reduction = c("cca", "rpca")[1])
sce <- IntegrateData(anchorset = datalist, dims = 1:40)
#ScaleData
DefaultAssay(sce) <- "integrated"

#sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000,
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
sce <- ScaleData(sce, features =  rownames(sce))
sce <- RunPCA(sce, features = VariableFeatures(sce)) 

dimplot <- DimPlot(sce, reduction = "pca",group.by = 'Samples') 
elbowplot <- ElbowPlot(sce, ndims=50, reduction="pca") 
sc_pca <- dimplot+elbowplot
sc_pca

Dims <- 30
#sce <- RunUMAP(sce, dims=1:Dims, reduction="pca")
###tsne 降维
sce <- RunTSNE(sce, 
               dims=1:Dims, 
               reduction="pca",
               perplexity=30,
               max_iter=1000)
umap.after<-DimPlot(sce,group.by='Samples',
                    reduction="tsne",
                    label = "T", 
                    pt.size = 0.2,
                    label.size = 0)+
  ggtitle('After batch elimination')+
  theme(legend.position = 'right')
umap.after
figs1ab<-mg_merge_plot(ump.raw,umap.after,nrow = 1,ncol = 2,labels = c('B','C'),widths = c(1,1.5))
figs1ab
head(summary_cells)
summary_cells$Samples=rownames(summary_cells)
summary_cells1=summary_cells
summary_cells1$fit=summary_cells1$raw_count-summary_cells$clean_count
summary_cells1=reshape2::melt(summary_cells1[,c("clean_count","Samples","fit")])
summary_cells1$variable=ifelse(summary_cells1$variable=='fit','Filtering number','Reserve number')
summary_cells1$variable=factor(summary_cells1$variable,levels = c('Filtering number','Reserve number'))
figs1c<-ggbarplot(summary_cells1, x = "Samples", y="value", color="black", fill="variable",
                  legend="right", 
                  legend.title="", main="Before and after filtering",
                  font.main = c(14,"bold", "black"), font.x = c(12, "bold"), 
                  font.y=c(12,"bold")) + 
  theme_bw() +
  rotate_x_text() + 
  labs(x = "", y = "Cell Number") + 
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title = element_text(face = "bold"), 
    plot.title = element_text(face = "bold"), 
    legend.title = element_text(face = "bold")) 
figs1abc<-mg_merge_plot(figs1c,figs1ab,labels = c('A',''),nrow = 2,ncol = 1,heights = c(1,1.5))
figs1abc
#聚类分析
DefaultAssay(sce) <- "integrated"
Resolution <- 0.9
sce <- FindNeighbors(object = sce, dims = 1:Dims)
sce <- FindClusters(object = sce, resolution = Resolution)
#先提取T细胞
DefaultAssay(sce) <- "RNA"
figs1d<-VlnPlot(sce,features = c('PTPRC','CD3D', 'CD3E', 'CD3G'),pt.size = 0,
                group.by = 'seurat_clusters',ncol = 2)+
  theme(legend.position = 'none')
figs1d
save(sce,file = 'sce.all.RData')
#重新聚类
load('sce.all.RData')
Idents(sce)='seurat_clusters'
sce<-subset(sce,idents =c(0,6))

DefaultAssay(sce) <- "integrated"

library(clustree)
sce <- FindNeighbors(sce, dims = 1:Dims)
sce <- FindClusters(
  object = sce,
  resolution = c(seq(.1,1,.1))
)
colnames(sce@meta.data)
clustree(sce@meta.data, prefix = "integrated_snn_res.")
table(sce$Samples)

pdf('results/ana/01.scRNA/clust.snn_res.pdf',he=15,wi=15)
clustree(sce@meta.data, prefix = "integrated_snn_res.")
dev.off()

save(sce,file = 'sce.T.RData')
load('sce.T.RData')
# #提取CD8A|CDB表达的
cd8_list <- list(c("CD8A","CD8B"))
DefaultAssay(sce) <- "RNA"
sce <- AddModuleScore(object = sce, features = cd8_list, name = "cd8")
FeaturePlot(object = sce, features = "cd81")
cell.name=rownames(sce@meta.data[which(sce@meta.data$cd81>0),])
cell.name
sce$cell.name=rownames(sce@meta.data)
Idents(sce)='cell.name'
sce.cd8<-subset(sce,idents =cell.name)
#sce.cd8=sce
table(sce.cd8$Samples)
DefaultAssay(sce.cd8) <- "integrated"
#cd8.marker=c('TNF','IL2','CXCR3','TBX21','IL4','IL5','CCR4','GATA3','IL9','IL10','IRF4','RORC','IL22')

Resolution <- 0.5
sce.cd8 <- FindNeighbors(object = sce.cd8, dims = 1:Dims)
sce.cd8 <- FindClusters(object = sce.cd8, resolution = Resolution)
DefaultAssay(sce.cd8) <- "RNA"
VlnPlot(sce.cd8,features = c('PTPRC','CD3D', 'CD3E', 'CD3G','CD8A','CD8B','CD4'),pt.size = 0,
        group.by = 'seurat_clusters',ncol = 4)+
  theme(legend.position = 'none')
figs1e<-VlnPlot(sce.cd8,features = c('PTPRC','CD3D', 'CD3E', 'CD3G','CD8A','CD8B'),pt.size = 0,
                group.by = 'seurat_clusters',ncol = 3)+
  theme(legend.position = 'none')
figs1e

#CD4和CD8和CD3
#CD4
#Th1:CXCR3、T-bet（TBX21）、TNFα（TNF）
#TH2:GATA-3(GATA3)、CD294(PTGDR2)、ST2
#Th17：CCR6、IL-17
#Th22:CCR10、IL-22、CCR6
#Tfh:CXCR5、IL-21、ICOS、PD-1
#Treg:CD25和CD127

#CD8
#Tc1（TNFα，INFγ，IL-2(IL2)， CXCR3，TBX21）；
#Tc2（IL-4，IL-5，CCR4，GATA3）；
#Tc9（IL-9，IL-10，IRF4）；
#Tc17（CCR6，KLRB1，IL-17，IRF4，RORC）
#Tc22（IL-2，IL-22，TNFα）
Idents(sce.cd8)='seurat_clusters'
#sce.cd8<-subset(sce.cd8,idents =c(0,1,2))
DefaultAssay(sce.cd8)<-'RNA'
VlnPlot(sce.cd8,features = c('CD8B','TNF','IL2','CXCR3','TBX21','IL4','IL5','CCR4','GATA3','IL9','IL10','IRF4','RORC','IL22'),pt.size = 0.1,
        group.by = 'seurat_clusters',ncol = 4)+
  theme(legend.position = 'none')

figs1de<-mg_merge_plot(figs1d,figs1e,nrow = 2,ncol = 1,labels = c('D','E'))
figs1de
figs1<-mg_merge_plot(figs1abc,figs1de,nrow = 1,ncol = 2,widths = c(2,1))
figs1
ggsave('results/ana/01.scRNA/FigS1.pdf',figs1,height = 10,width = 20)

#
sce.cd8$cell_type<-paste0('CD8_',sce.cd8$seurat_clusters)
table(sce.cd8$cell_type)
#marker基因的筛选
#寻找差异基因时的差异倍数
Logfc = 0.35
#差异基因时最小的表达比例
Minpct = 0.35
DefaultAssay(sce.cd8) <- "RNA"
Idents(sce.cd8)<-'cell_type'
sce.cd8.markers <- FindAllMarkers(object = sce.cd8,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
sce.cd8.markers["pct.diff"]=sce.cd8.markers$pct.1-sce.cd8.markers$pct.2
sce.cd8.markers <- sce.cd8.markers[sce.cd8.markers$p_val_adj<0.05,]
length(unique(sce.cd8.markers$gene))
head(sce.cd8.markers)
write.table(sce.cd8.markers,'results/files/CD8.scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')
write.table(sce.cd8.markers,'results/ana/01.scRNA/CD8.scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')

### 选择前5个marker基因
Top5 <- sce.cd8.markers %>% group_by(cluster) %>% slice_max(n =10, order_by = avg_logFC)  

Top5 <- intersect(unique(Top5$gene),rownames(sce.cd8@assays$RNA@meta.features))

sc_marker_dotplot <- DotPlot(object = sce.cd8, features = Top5,cols=c("blue", "red"),scale = T)+ 
  RotatedAxis()+ ggtitle("Top 10 Marker Genes")+ 
  theme(plot.title = element_text(hjust = 0.5)) +xlab('')+ylab('')

sc_marker_dotplot

fig1a<-DimPlot(sce.cd8,group.by = 'cell_type',
               reduction="tsne",
               label = "F", 
               pt.size = 1,
               label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
fig1a
as.matrix(table(sce.cd8$Samples,sce.cd8$cell_type))
length(unique(sce.cd8$Samples))
library("ggplot2")
cell_por<-as.matrix(table(sce.cd8$Samples,sce.cd8$cell_type))
cell_por=apply(cell_por,1,function(x){return(x/sum(x))})
cell_por=reshape2::melt(cell_por)
colnames(cell_por)<-c("cell_type",'Samples',"perc")

cell_por$proportion = paste(round(100*cell_por$perc ), "%")
cell_por=arrange(cell_por, perc) 
fig1b=ggplot(cell_por, aes(x = "", y = perc, fill = cell_type)) +
  geom_col(color = "black") +
  ggrepel::geom_label_repel(aes(label = proportion), color = c("white"),
                            position = position_stack(vjust = 0.5),
                            show.legend = FALSE) +
  guides(fill = guide_legend(title = "cell_type")) +
  scale_fill_hue() +
  coord_polar(theta = "y") + 
  theme_void()+facet_wrap(.~Samples,ncol = 2)+
  theme(strip.text=element_text(size = 10,family="Times",face="plain"))
fig1b
##cnv分析
# save(sce.cd8,file = 'sce.cd8.RData')
# load('sce.cd8.RData')
# library(copykat)
# copykat.test <- copykat(rawmat=sce.cd8@assays$RNA@counts,
#                         id.type="S",
#                         cell.line="no",
#                         ngene.chr=5,
#                         #每个染色体中至少有 5 个基因来计算 DNA 拷贝数
#                         win.size=25,
#                         #每个片段至少取 25 个基因
#                         KS.cut=0.15,
#                         #0-1,值越大灵敏度越低
#                         sam.name="LUAD",
#                         #随意固定一个名称
#                         distance="euclidean",
#                         n.cores=30
#                         #并行计算
# )
# save(copykat.test,file = 'copykat.test.RData')
# 
# #读取CNV的结果
# copykat.test<-read.delim('LUAD_copykat_prediction.txt',sep='\t',header = T)
# head(copykat.test)
# table(copykat.test$copykat.pred)
# copykat.test$copykat.pred=ifelse(copykat.test$copykat.pred=='aneuploid','malignant','no_malignant')
# table(copykat.test$copykat.pred)
# # malignant no_malignant 
# #155          154 
# rownames(copykat.test)=copykat.test$cell.names
# copykat.test=copykat.test[rownames(sce.cd8@meta.data),]
# #添加分组
# sce.cd8 <- AddMetaData(sce.cd8, copykat.test$copykat.pred,col.name = "copykat.pred")
# sce.cd8$copykat.pred[is.na(sce.cd8$copykat.pred)]<-'Unknown'
# table(sce.cd8$copykat.pred)
# # malignant no_malignant 
# # 155          154 
# as.matrix(table(sce.cd8$copykat.pred,sce.cd8$cell_type))
# library("ggplot2")
# cnv_por<-as.matrix(table(sce.cd8$cell_type,sce.cd8$copykat.pred))
# cnv_por=apply(cnv_por,1,function(x){return(x/sum(x))})
# cnv_por=reshape2::melt(cnv_por)
# colnames(cnv_por)<-c('cnv',"cluster","perc")
# 
# cnv_por$proportion = paste(round(100*cnv_por$perc ), "%")
# cnv_por=arrange(cnv_por, perc) 
# fig1c=ggplot(cnv_por, aes(x = "", y = perc, fill = cnv)) +
#   geom_col(color = "black") +
#   ggrepel::geom_label_repel(aes(label = proportion), color = c("white"),
#                             position = position_stack(vjust = 0.5),
#                             show.legend = FALSE) +
#   guides(fill = guide_legend(title = "cluster")) +
#   scale_fill_hue() +
#   coord_polar(theta = "y") + 
#   theme_void()+facet_wrap(.~cluster,ncol = 1)+
#   theme(strip.text=element_text(size = 10,family="Times",face="plain"))
# fig1c
#富集分析
#clust_name=as.character(unique(sce.markers$cluster))
library(clusterProfiler)
library(org.Hs.eg.db)
ids=bitr(sce.cd8.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers2=merge(sce.cd8.markers,ids,by.x='gene',by.y='SYMBOL')

## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers2$ENTREZID, sce.markers2$cluster) 
## KEGG
sce.markers2.enrich.res <- compareCluster(gcSample,
                                          fun = "enrichKEGG",
                                          organism = "hsa", pvalueCutoff = 0.05)
fig1d<-dotplot(sce.markers2.enrich.res)+ 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size=10),
        axis.text.y=element_text(size=10))
fig1d
fig1abc<-mg_merge_plot(fig1a,fig1b,nrow = 1,ncol = 2,labels = c('A','B'),widths = c(1,1))
fig1abc
fig1de<-mg_merge_plot(sc_marker_dotplot,fig1d,nrow = 1,ncol = 2,labels = c('C','D'),widths = c(1.5,1))

fig1<-mg_merge_plot(fig1abc,fig1de,nrow = 2,ncol = 1,heights = c(1,1),align = 'h')
ggsave('results/ana/01.scRNA/Fig1.pdf',fig1,height = 10,width = 15)


table(sce.cd8.markers$cluster)
marker1=unique(sce.cd8.markers$gene)
length(marker1)
VlnPlot(sce.cd8,features = marker1[91:100],pt.size = 0,group.by = 'cell_type')

cd8_0.marker1=c("ARRDC3","IFNG","FCGR3A","KLRF1","GZMA","FGFBP2","PYHIN1","CCL3","APMAP","CCL4L2","TPST2","CMC1","SPON2","CD47","RAC2","MYLIP","PLAC8","TAP1")
cd8_1.marker1=c("LMNA","RP11-138A9.1","SLC7A5","SQSTM1","JUND","AC016831.7","DUSP4","PLP2","SKIL","RP11-138A9.2","MYADM","MRPS6","FKBP4","RALGAPA1")
figs2<-VlnPlot(sce.cd8,features = c(cd8_0.marker1,cd8_1.marker1),pt.size = 0,
               group.by = 'cell_type',ncol = 6)
ggsave('results/ana/01.scRNA/FigS2.pdf',figs2,height = 15,width = 15)



marker.uniq<-rbind.data.frame(data.frame(cell_type='CD8_0',gene=cd8_0.marker1),
                              data.frame(cell_type='CD8_1',gene=cd8_1.marker1))
table(marker.uniq$cell_type)
write.table(marker.uniq,'results/ana/01.scRNA/marker.uniq.txt',quote = F,row.names = F,sep='\t')
write.table(marker.uniq,'results/files/marker.uniq.txt',quote = F,row.names = F,sep='\t')

save.image('project_00.RData')
##第二部分#####################
options(stringsAsFactors = F,check.bounds = F)

#2、Bulk-Seq数据的准备####
dir.create('results/ana/02.data.pre')
#2.1 TCGA数据的准备####
tcga_cli<-read.delim('origin_datas/TCGA/Merge_LUAD_clinical.txt',sep='\t',header = T,check.names = F)
tcga_cli=tcga_cli[,c("A0_Samples","A17_Age","A18_Sex","A3_T","A4_N","A5_M","A6_Stage","A1_OS","A2_Event","tobacco_smoking_history")]
colnames(tcga_cli)=c('Samples','Age','Gender','T.Stage','N.Stage','M.Stage','Stage','OS.time','OS','Smoking')
rownames(tcga_cli)=tcga_cli$Samples
tcga_cli=crbind2DataFrame(tcga_cli)
fivenum(tcga_cli$Age)
table(tcga_cli$Gender)

table(tcga_cli$T.Stage)
tcga_cli$T.Stage=gsub('[abc]','',tcga_cli$T.Stage)
tcga_cli$T.Stage[tcga_cli$T.Stage==''|tcga_cli$T.Stage=='TX']<-NA

table(tcga_cli$N.Stage)
tcga_cli$N.Stage[tcga_cli$N.Stage==''|tcga_cli$N.Stage=='NX']<-NA

table(tcga_cli$M.Stage)
tcga_cli$M.Stage=gsub('[abc]','',tcga_cli$M.Stage)
tcga_cli$M.Stage[tcga_cli$M.Stage==''|tcga_cli$M.Stage=='MX']<-NA

table(tcga_cli$Stage)
tcga_cli$Stage=gsub('Stage ','',tcga_cli$Stage)
tcga_cli$Stage=gsub('[ABC]','',tcga_cli$Stage)
tcga_cli$Stage[tcga_cli$Stage=='']<-NA


table(tcga_cli$OS)
tcga_cli$OS=ifelse(tcga_cli$OS=='Alive',0,ifelse(tcga_cli$OS=='Dead',1,NA))

tcga_cli$Samples=paste0(tcga_cli$Samples,'-01')
rownames(tcga_cli)=tcga_cli$Samples

table(tcga_cli$Smoking)
tcga_cli$OS.time
tcga_cli=tcga_cli[which(tcga_cli$OS.time>0),]
#去除10年以上
tcga_cli=tcga_cli[which(tcga_cli$OS.time<10 * 365),]


#表达谱
tcga_tmp<-read.delim('origin_datas/TCGA/Merge_TCGA-LUAD_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
rownames(tcga_tmp)[1:6]
#肿瘤样本和正常样本
sample_N<-colnames(tcga_tmp)[substr(colnames(tcga_tmp),14,15)=='11']
sample_T<-colnames(tcga_tmp)[substr(colnames(tcga_tmp),14,15)=='01']
sample_T=intersect(sample_T,tcga_cli$Samples)
length(sample_T)
#491
length(sample_N)
#59
table(gene.type$TYPE)
gene.type=gene.type[which(gene.type$TYPE=='protein_coding'),'SYMBOL']
tcga_tmp=tcga_tmp[intersect(rownames(tcga_tmp),gene.type),]


tcga_tmp_log2<-log2(tcga_tmp[,c(sample_N,sample_T)]+1)
tcga_tmp_log2_T<-tcga_tmp_log2[,sample_T]
tcga_cli=tcga_cli[sample_T,]
tcga.group=rbind.data.frame(data.frame(Samples=sample_T,Type='Tumor'),
                            data.frame(Samples=sample_N,Type='Normal'))

write.table(tcga.group,'results/ana/02.data.pre/tcga.group.txt',quote = F,row.names = F,sep='\t')
write.table(tcga_tmp_log2,'results/ana/02.data.pre/tcga_tpm_log.txt',quote = F,row.names = T,sep='\t')
write.table(tcga_cli,'results/ana/02.data.pre/tcga_cli.txt',quote = F,row.names = F,sep='\t')
write.table(tcga_tmp_log2_T,'results/ana/02.data.pre/tcga_tmp_log2_T.txt',quote = F,row.names = T,sep='\t')

dim(tcga_tmp_log2)
#25483   550
#2.2 GSE31210 ####
load('origin_datas/GEO/luad_data.RData')
GSE31210_cli<-luad_data$GSE31210$GSE31210_cli
head(GSE31210_cli)
GSE31210.group=GSE31210_cli[,c("Acc","tissue")]
table(GSE31210.group$tissue)
colnames(GSE31210.group)=c('Samples','Type')
GSE31210.group$Type=ifelse(GSE31210.group$Type=='normal lung','N','T')
table(GSE31210.group$Type)
# N   T 
# 20 226 
GSE31210_cli=GSE31210_cli[,c("Acc","death","days before death/censor","age (years)","gender","smoking status","pathological stage","gene alteration status")]
colnames(GSE31210_cli)=c('Samples','OS','OS.time','Age','Gender','Smoking','Stage','gene_muti')
GSE31210_cli$OS=gsub('dead',1,GSE31210_cli$OS)
GSE31210_cli$OS=gsub('alive',0,GSE31210_cli$OS)
head(GSE31210_cli)
rownames(GSE31210_cli)=GSE31210_cli$Samples
GSE31210_cli=GSE31210_cli[GSE31210.group[GSE31210.group$Type=='T','Samples'],]

GSE31210_exp<-luad_data$GSE31210$GSE31210_exp
GSE31210_exp[1:4,1:4]
GSE31210_exp=GSE31210_exp[GSE31210_exp$gene!='',-1]
GSE31210_exp[1:4,1:4]
range(GSE31210_exp)
GSE31210_exp=log2(GSE31210_exp+1)
GSE31210_exp=GSE31210_exp[,GSE31210_cli$Samples]
dim(GSE31210_exp)
#21655   226

range(GSE31210_exp)

table(GSE31210_cli$Gender)
GSE31210_cli$Gender=ifelse(GSE31210_cli$Gender=='female','Female','Male')

table(GSE31210_cli$Smoking)

table(GSE31210_cli$Stage)
GSE31210_cli$Stage=gsub('[AB]','',GSE31210_cli$Stage)

table(GSE31210_cli$gene_muti)
GSE31210_cli=crbind2DataFrame(GSE31210_cli)
write.table(GSE31210_exp,'results/ana/02.data.pre/GSE31210_exp.txt',quote = F,row.names = T,sep='\t')
write.table(GSE31210_cli,'results/ana/02.data.pre/GSE31210_cli.txt',quote = F,row.names = F,sep='\t')
#2.3 GSE50081#####
GSE50081_cli<-luad_data$GSE50081$GSE50081_cli
table(GSE50081_cli$histology)
GSE50081_cli=GSE50081_cli[which(GSE50081_cli$histology=='adenocarcinoma'),]
GSE50081_cli=data.frame(Samples=GSE50081_cli$Acc,
                        T.Stage=GSE50081_cli$`t-stage`,
                        N.Stage=GSE50081_cli$`n-stage`,
                        M.Stage=GSE50081_cli$`m-stage`,
                        Stage=GSE50081_cli$Stage,
                        Age=GSE50081_cli$age,
                        Gender=GSE50081_cli$Sex,
                        OS.time=GSE50081_cli$`survival time`,
                        OS=GSE50081_cli$status)
head(GSE50081_cli)
rownames(GSE50081_cli)=GSE50081_cli$Samples
table(GSE50081_cli$T.Stage)
GSE50081_cli$T.Stage=paste0('T',GSE50081_cli$T.Stage)
table(GSE50081_cli$N.Stage)
GSE50081_cli$N.Stage=paste0('N',GSE50081_cli$N.Stage)
table(GSE50081_cli$M.Stage)
GSE50081_cli$M.Stage=paste0('M',GSE50081_cli$M.Stage)
table(GSE50081_cli$Stage)
GSE50081_cli$Stage=gsub('[AB]','',GSE50081_cli$Stage)
GSE50081_cli$Stage=ifelse(GSE50081_cli$Stage==1,'I','II')
table(GSE50081_cli$Gender)
GSE50081_cli$Gender=ifelse(GSE50081_cli$Gender=='F','Female','Male')
table(GSE50081_cli$OS)
GSE50081_cli$OS=ifelse(GSE50081_cli$OS=='alive',0,1)
GSE50081_cli$Age
GSE50081_cli$OS.time
GSE50081_cli$OS.time=ceiling(GSE50081_cli$OS.time*365)
GSE50081_exp<-luad_data$GSE50081$GSE50081_exp
GSE50081_exp[1:4,1:4]
GSE50081_exp=GSE50081_exp[,-1]
range(GSE50081_exp)
GSE50081_cli=GSE50081_cli[which(GSE50081_cli$OS.time<10 * 365),]

GSE50081_exp=GSE50081_exp[,GSE50081_cli$Samples]
write.table(GSE50081_exp,'results/ana/02.data.pre/GSE50081_exp.txt',quote = F,row.names = T,sep='\t')
write.table(GSE50081_cli,'results/ana/02.data.pre/GSE50081_cli.txt',quote = F,row.names = F,sep='\t')
#2.4 GSE30219####
GSE30219_cli<-luad_data$GSE30219$GSE30219_cli

GSE30219_cli=data.frame(Samples=GSE30219_cli$Acc,
                        T.Stage=GSE30219_cli$`pt stage`,
                        N.Stage=GSE30219_cli$`pn stage`,
                        M.Stage=GSE30219_cli$`pm stage`,
                        Age=GSE30219_cli$`age at surgery`,
                        Gender=GSE30219_cli$gender,
                        OS.time=GSE30219_cli$`follow-up time (months)`,
                        OS=GSE30219_cli$status)
head(GSE30219_cli)
rownames(GSE30219_cli)=GSE30219_cli$Samples
table(GSE30219_cli$T.Stage)
GSE30219_cli$T.Stage[GSE30219_cli$T.Stage=='NTL'|GSE30219_cli$T.Stage=='TX']<-NA
table(GSE30219_cli$N.Stage)
GSE30219_cli$N.Stage[GSE30219_cli$N.Stage=='NTL'|GSE30219_cli$N.Stage=='NX']<-NA
table(GSE30219_cli$M.Stage)
GSE30219_cli$M.Stage[GSE30219_cli$M.Stage=='NTL'|GSE30219_cli$M.Stage=='MX']<-NA
table(GSE30219_cli$Gender)
GSE30219_cli$Gender=ifelse(GSE30219_cli$Gender=='F','Female',ifelse(GSE30219_cli$Gender=='M','Male',NA))
table(GSE30219_cli$OS)
GSE30219_cli$OS=ifelse(GSE30219_cli$OS=='ALIVE',0,ifelse(GSE30219_cli$OS=='DEAD',1,NA))
GSE30219_cli$OS.time
#去掉NA
library(tidyr)
GSE30219_cli <- GSE30219_cli %>% drop_na(OS.time)
GSE30219_cli=GSE30219_cli[GSE30219_cli$OS.time>0,]
GSE30219_cli$OS.time=ceiling(GSE30219_cli$OS.time*30)

GSE30219_exp<-luad_data$GSE30219$GSE30219_exp
GSE30219_exp[1:4,1:4]
GSE30219_exp=GSE30219_exp[,-1]
range(GSE30219_exp)
GSE30219_cli=GSE30219_cli[which(GSE30219_cli$OS.time<10 * 365),]

GSE30219_exp=GSE30219_exp[,GSE30219_cli$Samples]
write.table(GSE30219_exp,'results/ana/02.data.pre/GSE30219_exp.txt',quote = F,row.names = T,sep='\t')
write.table(GSE30219_cli,'results/ana/02.data.pre/GSE30219_cli.txt',quote = F,row.names = F,sep='\t')
#2.5 GSE19188 #####
GSE19188_cli<-luad_data$GSE19188$GSE19188_cli
GSE19188_cli=data.frame(Samples=GSE19188_cli$Acc,
                        OS.time=GSE19188_cli$`overall survival`,
                        OS=GSE19188_cli$status,
                        Gender=GSE19188_cli$gender)
head(GSE19188_cli)
rownames(GSE19188_cli)=GSE19188_cli$Samples
table(GSE19188_cli$OS)
GSE19188_cli=GSE19188_cli[which(GSE19188_cli$OS != 'Not available'),]
GSE19188_cli$OS=ifelse(GSE19188_cli$OS=='alive',0,1)
table(GSE19188_cli$Gender)
GSE19188_cli$Gender=ifelse(GSE19188_cli$Gender=='F','Female','Male')
GSE19188_cli$OS.time
GSE19188_cli$OS.time=ceiling(as.numeric(GSE19188_cli$OS.time)*30)

GSE19188_exp<-luad_data$GSE19188$GSE19188_exp
GSE19188_exp[1:4,1:4]
GSE19188_exp=GSE19188_exp[,-1]
range(GSE19188_exp)
GSE19188_exp=log2(2^GSE19188_exp+1)
range(GSE19188_exp)

GSE19188_cli=GSE19188_cli[which(GSE19188_cli$OS.time<10 * 365),]
GSE19188_exp=GSE19188_exp[,GSE19188_cli$Samples]
write.table(GSE19188_exp,'results/ana/02.data.pre/GSE19188_exp.txt',quote = F,row.names = T,sep='\t')
write.table(GSE19188_cli,'results/ana/02.data.pre/GSE19188_cli.txt',quote = F,row.names = F,sep='\t')
#2.6 GSE37745 #### 
GSE37745_cli<-luad_data$GSE37745$GSE37745_cli
GSE37745_cli=data.frame(Samples=GSE37745_cli$Acc,
                        Stage=GSE37745_cli$`tumor stage`,
                        OS=GSE37745_cli$dead,
                        OS.time=GSE37745_cli$`days to determined death status`)
head(GSE37745_cli)
table(GSE37745_cli$OS)
GSE37745_cli$OS=ifelse(GSE37745_cli$OS=='yes',1,0)
GSE37745_cli$OS.time
table(GSE37745_cli$Stage)
GSE37745_cli$Stage=gsub('[ab]','',GSE37745_cli$Stage)
GSE37745_cli$Stage=ifelse(GSE37745_cli$Stage==1,'I',ifelse(GSE37745_cli$Stage==2,'II',ifelse(GSE37745_cli$Stage==3,'III','IV')))
GSE37745_exp<-luad_data$GSE37745$GSE37745_exp
GSE37745_exp[1:4,1:5]
GSE37745_exp=GSE37745_exp[,-1]
range(GSE37745_exp)

GSE37745_cli=GSE37745_cli[which(GSE37745_cli$OS.time<10 * 365),]
GSE37745_exp=GSE37745_exp[,GSE37745_cli$Samples]
write.table(GSE37745_exp,'results/ana/02.data.pre/GSE37745_exp.txt',quote = F,row.names = T,sep='\t')
write.table(GSE37745_cli,'results/ana/02.data.pre/GSE37745_cli.txt',quote = F,row.names = F,sep='\t')
#3、CD8 T细胞亚群的进一步筛选#########
dir.create('results/ana/03.key.cell')
marker.uniq<-read.delim('results/ana/01.scRNA/marker.uniq.txt',sep='\t',header = T)
ssGSEAScore_by_genes=function(gene.exp,genes){
  #library('GSVA')
  #library(GSEABase)
  #all.list=list()
  gs=GSEABase::GeneSet(setName='GeneSet', setIdentifier=paste0("101"),geneIds=unique(genes),GSEABase::SymbolIdentifier()) 
  
  gsc <- GSEABase::GeneSetCollection(list(gs))
  fl <- tempfile()
  GSEABase::toGmt(gsc, fl)
  cgeneset=GSEABase::getGmt(fl)
  ssGSEA.geneset <- GSVA::gsva(as.matrix(gene.exp), cgeneset,method='ssgsea',
                               min.sz=1, max.sz=5000, verbose=TRUE)
  #detach('package:GSVA')
  #detach('package:GSEABase')
  #row.names(ssGSEA.geneset)
  return(ssGSEA.geneset)
}
cell.score<-function(exp,gene){
  pathway_score<-data.frame()
  for (i in unique(gene[,1])){
    print(i)
    gene_set=gene[gene[,1]==i,"gene"]
    score=ssGSEAScore_by_genes(exp,gene_set)
    rownames(score)=i
    pathway_score=rbind.data.frame(pathway_score,score)
  }
  return(t(pathway_score))
}
table(marker.uniq$cell_type)
marker.uniq=rbind.data.frame(marker.uniq,
                             data.frame(cell_type='CD8_0',gene=c('PTPRC','CD3D', 'CD3E', 'CD3G','CD8A','CD8B')),
                             data.frame(cell_type='CD8_1',gene=c('PTPRC','CD3D', 'CD3E', 'CD3G','CD8A','CD8B')))
tcga.score<-cell.score(exp = tcga_tmp_log2,gene = marker.uniq)
head(tcga.score)
fig2a<-mg_PlotMutiBoxplot(data = tcga.score[tcga.group$Samples,],
                          group = tcga.group$Type,
                          group_cols = c('blue','red'),
                          test_method = 'wilcox.test',
                          xangle = 30,
                          add = 'boxplot', 
                          legend.pos = 'top',ylab = 'CD8 T cell Score')
fig2a
#森林图
tcga.cell.cox <- cox_batch(t(tcga.score[tcga_cli$Samples,]),
                           time = tcga_cli$OS.time,
                           event = tcga_cli$OS)
tcga.cell.cox
# tcga.cell.cox$Name=rownames(tcga.cell.cox)
# tcga.cell.cox <- data.frame(Names=rownames(tcga.cell.cox),
#                             p.value=tcga.cell.cox$p.value,
#                             tcga.cell.cox$HR,
#                             tcga.cell.cox$`Low 95%CI`,
#                             tcga.cell.cox$`High 95%CI`)
# pdf('results/ana/03.key.cell/fig2b.pdf',height = 5,width = 7,onefile = F)
# mg_forestplot_v2(tcga.cell.cox,xlog = T,colgap = 8,lineheight = 10)
# dev.off()
#KM曲线
library(survminer)
biosurvival<-function(OS,OS.time,score,palette = ggsci::pal_lancet()(9)[c(2,3)],title){
  dat=data.frame(OS.time=OS.time,
                 OS=OS,
                 score=as.numeric(score))
  res.cut <- surv_cutpoint(dat,
                           time = "OS.time", 
                           event = "OS", 
                           variables = c("score"))
  
  cut_va=as.numeric(res.cut$cutpoint[1])
  p<-ggplotKMCox(dat = data.frame(dat$OS.time,
                                  dat$OS,
                                      ifelse(dat$score>median(dat$score),'High','low')),
                                      #ifelse(dat$score>cut_va,'High','low')),
                     add='',palette = palette,
                     show_confint = T,
                 labs = c('High','Low'),title = title)
  return(p)
}

fig2b<-biosurvival(OS =tcga_cli$OS,
                   OS.time = tcga_cli$OS.time/365,
                   score = as.numeric(tcga.score[tcga_cli$Samples,"CD8_0"]),
                   palette = ggsci::pal_lancet()(9)[c(1,3)],
                   title = 'CD8_0 Group')
fig2b
fig2c<-biosurvival(OS =tcga_cli$OS,
                   OS.time = tcga_cli$OS.time/365,
                   score = as.numeric(tcga.score[tcga_cli$Samples,"CD8_1"]),
                   palette = ggsci::pal_lancet()(9)[c(1,3)],
                   title = 'CD8_1 Group')
fig2c


#ROC曲线
# fig2d<-ggplotTimeROC(time = tcga_cli$OS.time/365,status = tcga_cli$OS,
#                      score = as.numeric(tcga.score[tcga_cli$Samples,"CD8_2"]),
#                      mks = c(1,2,3,4,5))+
#   ggsci::scale_color_lancet()
# fig2d
#TIDE
tcga_tide_dat <- t(scale(t(tcga_tmp_log2_T),scale = F))
write.table(tcga_tide_dat,
            file = 'results/ana/03.key.cell/tcga_tide_dat.txt',
            quote = F, sep = '\t')

tcga.tide<-read.csv('results/ana/03.key.cell/tcga.tide.res.csv',row.names = 1,stringsAsFactors = F)
head(tcga.tide)
tcga.tide=crbind2DataFrame(tcga.tide)
tcga.tide.cell=cbind.data.frame(tcga.tide,
                                tcga.score[rownames(tcga.tide),])
head(tcga.tide.cell)
sig_boxplot<-function(dat,leg,ylab,palette=ggsci::pal_lancet()(10)[3:4]){
  library(ggpubr)
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggboxplot(dat, 
               x='group', y='gene', color = 'group',
               palette =  palette,
               short.panel.labs = T,outlier.shape = NA)+
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab('')+labs(color=leg)
  return(pp)
}
tcga.score1<-cbind.data.frame(tcga.score[rownames(tcga.tide),],
                              TIDE=tcga.tide$TIDE)
# fig2d<-cor_point(x=tcga.score1$CD8_0,y = tcga.score1$TIDE,method = 'spearman')
# fig2d
# 
# fig2e<-cor_point(x=tcga.score1$CD8_1,y = tcga.score1$TIDE,method = 'spearman')
# fig2e

# mg_PlotMutiBoxplot(data=tcga.score[rownames(tcga.tide),],
#                           group = tcga.tide$Responder,
#                           group_cols =  ggsci::pal_lancet()(9)[c(1,3)],
#                           test_method = 'wilcox.test',
#                           xangle = 30,
#                           add = 'boxplot', 
#                           legend.pos = 'top',
#                           ylab = 'CD8 T cell Score')
# fig2e

fig2<-mg_merge_plot(fig2a,fig2b,fig2c,nrow = 1,ncol = 3,labels = c('A','B','C'))

ggsave('results/ana/03.key.cell/Fig2.pdf',fig2,height = 5,width = 15)

#4、WGCNA筛选关键的基因#####
dir.create('results/ana/04.WGCNA')
library(WGCNA)
?allowWGCNAThreads
allowWGCNAThreads(nThreads = 36)#允许R语言程序最大线程运行
enableWGCNAThreads(nThreads = 36)# 打开多线程
wgcna.cli <- tcga.score[,"CD8_0",drop=F]
#过滤
range(tcga_tmp_log2_T)
#wgcna.input=tcga_tmp_log2_T[apply(tcga_tmp_log2_T,1,function(x){return(sd(x)>0.5)}),]
wgcna.input <- tcga_tmp_log2_T[which(apply(tcga_tmp_log2_T,1,function(x){return(sum(x>=2.5))})>0.5*ncol(tcga_tmp_log2_T)),]
dim(wgcna.input)

datExpr=as.data.frame(t(wgcna.input))
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
dim(datExpr)

sampleTree = hclust(dist(datExpr), method = "complete")
pdf('results/ana/04.WGCNA/Fig3a.pdf',width = 6,height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")
# abline(h=140,col='red')
dev.off()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, 
                        powerVector = powers, 
                        verbose = 3)


pdf('results/ana/04.WGCNA/Fig3b.pdf',width = 8,height = 6)
cex1 = 0.85
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
softPower=sft$powerEstimate
softPower#10

net = blockwiseModules(datExpr, power = softPower, maxBlockSize = 22000,
                       TOMType = "unsigned", minModuleSize = 70,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       # saveTOMs = TRUE,
                       # saveTOMFileBase = "TPM-TOM-40", deepSplit = 3,
                       verbose = 3)
save(net, file = 'wgcna-net.RData')
table(net$colors)
mergedColors = labels2colors(net$colors)
table(mergedColors)
length(table(mergedColors))
pdf('results/ana/04.WGCNA/Fig3c.pdf',width = 8,height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    groupLabels = c("Module colors","GS.weight"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
wgcna.cli=crbind2DataFrame(wgcna.cli)
MEs_col = net$MEs[rownames(wgcna.cli), ]
colnames(MEs_col) = labels2colors(as.numeric(gsub('ME','',colnames(net$MEs))))

spms=wgcna.cli

dim(spms)[1]
modTraitCor = cor(MEs_col, as.data.frame(spms), use = "p")
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf('results/ana/04.WGCNA/Fig3d.pdf',width = 4,height = 8)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(spms), 
               yLabels = colnames(MEs_col), 
               cex.lab = 1, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
MEs = moduleEigengenes(datExpr[rownames(wgcna.cli),], mergedColors)$eigengenes
MEs = orderMEs(MEs)
colnames(MEs)=gsub('^ME','',colnames(MEs))
me.inds=which(colnames(MEs)!='grey')
if(length(me.inds)==1){
  cns=colnames(MEs)[me.inds]
  MEs=as.data.frame(MEs[,me.inds])
  colnames(MEs)=cns
}else if(length(me.inds)>1){
  MEs=MEs[,me.inds]
}


nGenes = ncol(datExpr)
nSamples = nrow(datExpr[rownames(wgcna.cli),])
geneTraitSignificance = as.data.frame(cor(datExpr[rownames(wgcna.cli),],spms, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(MEs)));
names(geneTraitSignificance) = paste("GS.", colnames(spms), sep="");
names(GSPvalue) = paste("p.GS.", colnames(spms), sep="")
geneModuleMembership = as.data.frame(cor(datExpr[rownames(wgcna.cli),], MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", names(geneModuleMembership), sep="");
names(MMPvalue) = paste("p.MM", names(geneModuleMembership), sep="");
write.table(cbind(Genes=row.names(geneTraitSignificance),geneTraitSignificance,GSPvalue,Module=mergedColors)
            ,file = "results/files/GeneTraitSignificance.txt",sep = '\t',quote = F,row.names = F)


# brown   模块
module_genes <- colnames(datExpr)[which(mergedColors=='brown')]
length(module_genes)#710

# 功能富集分析
com_genes_res <- enrichmentORA(module_genes,
                               mp_dbs=c('pathway_KEGG',
                                        'geneontology_Biological_Process',
                                        'geneontology_Cellular_Component',
                                        'geneontology_Molecular_Function'))
com_genes_res_filtered <- com_genes_res[com_genes_res$FDR < 0.05, ]
table(com_genes_res_filtered$DB)
# geneontology_Biological_Process geneontology_Cellular_Component geneontology_Molecular_Function 
#  858                              88                              85 
# pathway_KEGG 
# 58
pdf('results/ana/04.WGCNA/Fig4.pdf', width = 10, height = 8)
dotplot_batch(com_genes_res_filtered,
              dbs =c('geneontology_Biological_Process',
                     'geneontology_Cellular_Component',
                     'geneontology_Molecular_Function',
                     'pathway_KEGG'),top=10, FDR = T)
dev.off()
#hub基因的筛选
geneTraitSignificance1=geneTraitSignificance[module_genes,,drop=F]
geneModuleMembership1=geneModuleMembership[module_genes,,drop=F]
dim(geneTraitSignificance1[abs(geneTraitSignificance1$GS.CD8_2) >0.3,,drop=F])
dim(geneModuleMembership1[geneModuleMembership1$MMred>0.6,])

hub_gene=intersect(rownames(geneTraitSignificance1[abs(geneTraitSignificance1$GS.CD8_0)>0.6,,drop=F]),
                   rownames(geneModuleMembership1[geneModuleMembership1$MMbrown>0.8,,drop=F]))
length(hub_gene)#106
#绘图
pdf('results/ana/04.WGCNA/Fig3e.pdf',he=7,wi=7)
verboseScatterplot(abs(geneModuleMembership1[module_genes, 'MMbrown']),
                   abs(geneTraitSignificance1[module_genes, 'GS.CD8_0']),
                   xlab = paste("Module Membership in brown module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'brown')
abline(h=0.6,v=0.8,col='red')
dev.off()
#
write.table(hub_gene,file = 'results/ana/04.WGCNA/hub_gene.txt',quote = F,sep='\t',row.names = F,col.names = F)
#单因素cox分析
tcga.hub.cox <- cox_batch(tcga_tmp_log2_T[hub_gene,tcga_cli$Samples],
                          time = tcga_cli$OS.time,
                          event = tcga_cli$OS)
tcga.hub.cox
tcga.hub.cox$Name=rownames(tcga.hub.cox)
cluster.sig.gene=tcga.hub.cox[which(tcga.hub.cox$p.value<0.05),"Name"]
length(cluster.sig.gene)#42
write.table(tcga.hub.cox,'results/ana/04.WGCNA/tcga.hub.cox.txt',quote = F,row.names = F,sep='\t')
write.table(tcga.hub.cox,'results/files/tcga.hub.cox.txt',quote = F,row.names = F,sep='\t')

#5、一致性聚类####
dir.create('results/ana/05.subtype')
library(ConsensusClusterPlus)

clusterAlg_use=c('pam','hc','km','kmdist')[1]
distance_use=c('pearson','spearman','euclidean','canberra')[3]
#TCGA
df_exp=as.matrix(tcga_tmp_log2_T[cluster.sig.gene,])
#df_exp=sweep(df_exp,1,apply(df_exp, 1, median))
#df_exp=t(scale(t(df_exp)))
#df_exp=t(scale(t(df_exp),scale = F))
dim(df_exp)
clust_subtype_TCGA = ConsensusClusterPlus(df_exp
                                          , maxK = 10, reps = 500, pItem = 0.8
                                          , pFeature =1
                                          , title = "TCGA_subtype"
                                          , clusterAlg = clusterAlg_use
                                          , distance = distance_use
                                          , innerLinkage = 'ward.D2'
                                          , plot = "pdf"
                                          , seed = 123456
                                          , writeTable = T)

k=3
tcga.subtype=data.frame(Cluster=clust_subtype_TCGA[[k]]$consensusClass)
rownames(tcga.subtype)<-colnames(df_exp)
tcga.subtype$Cluster=paste0('MC',tcga.subtype$Cluster)
tcga.subtype[which(tcga.subtype$Cluster=='MC1'),1]='TMC3'
tcga.subtype[which(tcga.subtype$Cluster=='MC2'),1]='TMC2'
tcga.subtype[which(tcga.subtype$Cluster=='MC3'),1]='TMC1'

tcga.subtype$Cluster=gsub('TMC','clust',tcga.subtype$Cluster)
table(tcga.subtype$Cluster)
tcga.subtype$Samples=rownames(tcga.subtype)

ggplotKMCox(data.frame(time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                       , event = tcga_cli[rownames(tcga.subtype),]$OS
                       , tcga.subtype$Cluster)
            , add_text = ''
            #,labs = c('clust1','clust2')
)
# #GSE31210
# df_exp=as.matrix(GSE31210_exp[intersect(cluster.sig.gene,rownames(GSE31210_exp)),])
# #df_exp=sweep(df_exp,1,apply(df_exp, 1, median))
# #df_exp=t(scale(t(df_exp)))
# df_exp=t(scale(t(df_exp),scale = F))
# dim(df_exp)
# clust_subtype_gse31210 = ConsensusClusterPlus(df_exp
#                                           , maxK = 10, reps = 500, pItem = 0.8
#                                           , pFeature =1
#                                           , title = "GSE31210_subtype"
#                                           , clusterAlg = clusterAlg_use
#                                           , distance = distance_use
#                                           , innerLinkage = 'ward.D2'
#                                           , plot = "pdf"
#                                           , seed = 123456
#                                           , writeTable = T)
# 
# k=3
# gse31210.subtype=data.frame(Cluster=clust_subtype_gse31210[[k]]$consensusClass)
# rownames(gse31210.subtype)<-colnames(df_exp)
# gse31210.subtype$Cluster=paste0('MC',gse31210.subtype$Cluster)
# gse31210.subtype[which(gse31210.subtype$Cluster=='MC1'),1]='TMC1'
# gse31210.subtype[which(gse31210.subtype$Cluster=='MC2'),1]='TMC2'
# gse31210.subtype[which(gse31210.subtype$Cluster=='MC3'),1]='TMC3'
# 
# gse31210.subtype$Cluster=gsub('TMC','clust',gse31210.subtype$Cluster)
# table(gse31210.subtype$Cluster)
# gse31210.subtype$Samples=rownames(gse31210.subtype)
# 
# ggplotKMCox(data.frame(time = GSE31210_cli[rownames(gse31210.subtype),]$OS.time/365
#                        , event = GSE31210_cli[rownames(gse31210.subtype),]$OS
#                        , gse31210.subtype$Cluster)
#             , add_text = ''
#             #,labs = c('clust1','clust2')
# )
# #GSE19188
df_exp=as.matrix(GSE19188_exp[intersect(cluster.sig.gene,rownames(GSE19188_exp)),])
#df_exp=sweep(df_exp,1,apply(df_exp, 1, median))
df_exp=t(scale(t(df_exp)))
#df_exp=t(scale(t(df_exp),scale = F))
dim(df_exp)
clust_subtype_gse19188 = ConsensusClusterPlus(df_exp
                                              , maxK = 10, reps = 500, pItem = 0.8
                                              , pFeature =1
                                              , title = "GSE19188_subtype"
                                              , clusterAlg = clusterAlg_use
                                              , distance = distance_use
                                              , innerLinkage = 'ward.D2'
                                              , plot = "pdf"
                                              , seed = 123456
                                              , writeTable = T)

k=3
gse19188.subtype=data.frame(Cluster=clust_subtype_gse19188[[k]]$consensusClass)
rownames(gse19188.subtype)<-colnames(df_exp)
gse19188.subtype$Cluster=paste0('MC',gse19188.subtype$Cluster)
gse19188.subtype[which(gse19188.subtype$Cluster=='MC1'),1]='TMC1'
gse19188.subtype[which(gse19188.subtype$Cluster=='MC2'),1]='TMC2'
gse19188.subtype[which(gse19188.subtype$Cluster=='MC3'),1]='TMC3'

gse19188.subtype$Cluster=gsub('TMC','clust',gse19188.subtype$Cluster)
table(gse19188.subtype$Cluster)
gse19188.subtype$Samples=rownames(gse19188.subtype)

ggplotKMCox(data.frame(time = GSE19188_cli[rownames(gse19188.subtype),]$OS.time/365
                       , event = GSE19188_cli[rownames(gse19188.subtype),]$OS
                       , gse19188.subtype$Cluster)
            , add_text = ''
            #,labs = c('clust1','clust2')
)
#GSE30219
df_exp=as.matrix(GSE30219_exp[intersect(cluster.sig.gene,rownames(GSE30219_exp)),])
#df_exp=sweep(df_exp,1,apply(df_exp, 1, median))
df_exp=t(scale(t(df_exp)))
#df_exp=t(scale(t(df_exp),scale = F))
dim(df_exp)
clust_subtype_gse30219 = ConsensusClusterPlus(df_exp
                                              , maxK = 10, reps = 500, pItem = 0.8
                                              , pFeature =1
                                              , title = "GSE30219_subtype"
                                              , clusterAlg = clusterAlg_use
                                              , distance = distance_use
                                              , innerLinkage = 'ward.D2'
                                              , plot = "pdf"
                                              , seed = 123456
                                              , writeTable = T)

k=3
gse30219.subtype=data.frame(Cluster=clust_subtype_gse30219[[k]]$consensusClass)
rownames(gse30219.subtype)<-colnames(df_exp)
gse30219.subtype$Cluster=paste0('MC',gse30219.subtype$Cluster)
gse30219.subtype[which(gse30219.subtype$Cluster=='MC1'),1]='TMC3'
gse30219.subtype[which(gse30219.subtype$Cluster=='MC2'),1]='TMC1'
gse30219.subtype[which(gse30219.subtype$Cluster=='MC3'),1]='TMC2'

gse30219.subtype$Cluster=gsub('TMC','clust',gse30219.subtype$Cluster)
table(gse30219.subtype$Cluster)
gse30219.subtype$Samples=rownames(gse30219.subtype)

ggplotKMCox(data.frame(time = GSE30219_cli[rownames(gse30219.subtype),]$OS.time/365
                       , event = GSE30219_cli[rownames(gse30219.subtype),]$OS
                       , gse30219.subtype$Cluster)
            , add_text = ''
            #,labs = c('clust1','clust2')
)
# #GSE37745
# df_exp=as.matrix(GSE37745_exp[intersect(cluster.sig.gene,rownames(GSE37745_exp)),])
# #df_exp=sweep(df_exp,1,apply(df_exp, 1, median))
# #df_exp=t(scale(t(df_exp)))
# df_exp=t(scale(t(df_exp),scale = F))
# dim(df_exp)
# clust_subtype_gse37745 = ConsensusClusterPlus(df_exp
#                                               , maxK = 10, reps = 500, pItem = 0.8
#                                               , pFeature =1
#                                               , title = "GSE37745_subtype"
#                                               , clusterAlg = clusterAlg_use
#                                               , distance = distance_use
#                                               , innerLinkage = 'ward.D2'
#                                               , plot = "pdf"
#                                               , seed = 123456
#                                               , writeTable = T)
# 
# k=3
# gse37745.subtype=data.frame(Cluster=clust_subtype_gse37745[[k]]$consensusClass)
# rownames(gse37745.subtype)<-colnames(df_exp)
# gse37745.subtype$Cluster=paste0('MC',gse37745.subtype$Cluster)
# gse37745.subtype[which(gse37745.subtype$Cluster=='MC1'),1]='TMC1'
# gse37745.subtype[which(gse37745.subtype$Cluster=='MC2'),1]='TMC2'
# gse37745.subtype[which(gse37745.subtype$Cluster=='MC3'),1]='TMC3'
# 
# gse37745.subtype$Cluster=gsub('TMC','clust',gse37745.subtype$Cluster)
# table(gse37745.subtype$Cluster)
# gse37745.subtype$Samples=rownames(gse37745.subtype)
# rownames(GSE37745_cli)=GSE37745_cli$Samples
# ggplotKMCox(data.frame(time = GSE37745_cli[rownames(gse37745.subtype),]$OS.time/365
#                        , event = GSE37745_cli[rownames(gse37745.subtype),]$OS
#                        , gse37745.subtype$Cluster)
#             , add_text = ''
#             #,labs = c('clust1','clust2')
# )
#GSE50081
# df_exp=as.matrix(GSE50081_exp[intersect(cluster.sig.gene,rownames(GSE50081_exp)),])
# #df_exp=sweep(df_exp,1,apply(df_exp, 1, median))
# #df_exp=t(scale(t(df_exp)))
# df_exp=t(scale(t(df_exp),scale = F))
# dim(df_exp)
# clust_subtype_gse50081 = ConsensusClusterPlus(df_exp
#                                               , maxK = 10, reps = 500, pItem = 0.8
#                                               , pFeature =1
#                                               , title = "GSE50081_subtype"
#                                               , clusterAlg = clusterAlg_use
#                                               , distance = distance_use
#                                               , innerLinkage = 'ward.D2'
#                                               , plot = "pdf"
#                                               , seed = 123456
#                                               , writeTable = T)
# 
# k=3
# gse50081.subtype=data.frame(Cluster=clust_subtype_gse50081[[k]]$consensusClass)
# rownames(gse50081.subtype)<-colnames(df_exp)
# gse50081.subtype$Cluster=paste0('MC',gse50081.subtype$Cluster)
# gse50081.subtype[which(gse50081.subtype$Cluster=='MC1'),1]='TMC1'
# gse50081.subtype[which(gse50081.subtype$Cluster=='MC2'),1]='TMC2'
# gse50081.subtype[which(gse50081.subtype$Cluster=='MC3'),1]='TMC3'
# 
# gse50081.subtype$Cluster=gsub('TMC','clust',gse50081.subtype$Cluster)
# table(gse50081.subtype$Cluster)
# gse50081.subtype$Samples=rownames(gse50081.subtype)
# 
# ggplotKMCox(data.frame(time = GSE50081_cli[rownames(gse50081.subtype),]$OS.time/365
#                        , event = GSE50081_cli[rownames(gse50081.subtype),]$OS
#                        , gse50081.subtype$Cluster)
#             , add_text = ''
#             #,labs = c('clust1','clust2')
# )


fig4d<-ggplotKMCox(data.frame(time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                       , event = tcga_cli[rownames(tcga.subtype),]$OS
                       , tcga.subtype$Cluster)
            , add_text = ''
            ,palette = ggsci::pal_nejm()(9)[2:4]
            ,labs = c('clust1','clust2','clust3')
            ,title = 'TCGA'
)
fig4d

fig4e<-ggplotKMCox(data.frame(time = GSE30219_cli[rownames(gse30219.subtype),]$OS.time/365
                       , event = GSE30219_cli[rownames(gse30219.subtype),]$OS
                       , gse30219.subtype$Cluster)
                   , add_text = ''
                   ,palette = ggsci::pal_nejm()(9)[2:4]
                   ,labs = c('clust1','clust2','clust3')
                   ,title = 'GSE30219'
)
fig4e
write.table(data.frame(Samples=rownames(tcga.subtype),
                       time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                       , event = tcga_cli[rownames(tcga.subtype),]$OS
                       , tcga.subtype$Cluster),
            'results/ana/05.subtype/tcga.subtype.txt',quote = F,row.names = F,sep='\t')
write.table(data.frame(Samples=rownames(tcga.subtype),
                       time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                       , event = tcga_cli[rownames(tcga.subtype),]$OS
                       , tcga.subtype$Cluster),
            'results/files//tcga.subtype.txt',quote = F,row.names = F,sep='\t')
#
write.table(data.frame(Samples=rownames(gse30219.subtype),
                       time = GSE30219_cli[rownames(gse30219.subtype),]$OS.time/365
                       , event = GSE30219_cli[rownames(gse30219.subtype),]$OS
                       , gse30219.subtype$Cluster),
            'results/ana/05.subtype/gse30219.subtype.txt',quote = F,row.names = F,sep='\t')
write.table(data.frame(Samples=rownames(gse30219.subtype),
                       time = GSE30219_cli[rownames(gse30219.subtype),]$OS.time/365
                       , event = GSE30219_cli[rownames(gse30219.subtype),]$OS
                       , gse30219.subtype$Cluster),
            'results/files/gse30219.subtype.txt',quote = F,row.names = F,sep='\t')

fig4de<-mg_merge_plot(fig4d,fig4e,ncol = 2,nrow = 1,labels = c('D','E'))
fig4de
ggsave('results/ana/05.subtype/Fig4de.pdf',fig4de,height = 7,width = 12)
#6、临床表型的差异####
dir.create('results/ana/06.subtype.cli')
tcga_subtype.cli=merge(tcga.subtype,tcga_cli,by='Samples')
rownames(tcga_subtype.cli)=tcga_subtype.cli$Samples
tcga_subtype.cli$Event=ifelse(tcga_subtype.cli$OS==0,'Alive','Dead')
fivenum(na.omit(tcga_subtype.cli$Age))
tcga_subtype.cli$Age1=ifelse(tcga_subtype.cli$Age>60,'>60',ifelse(tcga_subtype.cli$Age<=60,'<=60',NA))
write.table(tcga_subtype.cli,'results/ana/06.subtype.cli/tcga_subtype.cli.txt',quote = F,row.names = F,sep='\t')

pie_compare_plot<-function(dat,gname,group_cols){
  library(dplyr)
  g_n<-as.character(unique(dat[,gname]))
  vname <- setdiff(colnames(dat), gname)
  pie.gname=data.frame()
  fisher.p <- c()
  for (i in vname) {
    tmp <- table(dat[,gname], dat[,i])
    p <- round(chisq.test(tmp)$p.value,4)
    names(p) <- i
    fisher.p <- c(fisher.p, p)
    pie.dat <- 
      tmp %>% as.data.frame() %>% group_by(Var1) %>% mutate(Pct = Freq/sum(Freq)) %>% as.data.frame()
    pie.gname=rbind.data.frame(pie.gname,pie.dat)
  }
  #plot
  vname_col<-list()
  for(i in 1:length(vname)){
    col_i<-length(as.character(unique(na.omit(dat[,vname[i]]))))
    vname_col[[i]] <- alpha(group_cols[i], 
                            sample(x = 1:10,size = col_i,replace = F )/10)
  }
  
  #names(vname_col)=vname
  col_num<-ncol(dat)
  row_num<-length(g_n)
  row_nums=1+2*row_num+2
  #c(1:col_num,2*row_num*col_num+1)
  #第一行
  nums1=c()
  for (i in 1:col_num){
    nums1=c(nums1,rep(i,3))
  }
  #两行
  nums2=c()
  for (j in 1:row_num){
    nums21=c()
    for (i in (col_num*j+1):((1+j)*col_num)){
      nums21=c(nums21,rep(i,3))
    }
    nums2=c(nums2,nums21,nums21)
  }
  
  #倒数第二行
  nums3=c()
  #(1+row_num)*col_num+1,(2+row_num)*col_num
  for (i in (((1+row_num)*col_num+1):((2+row_num)*col_num))){
    nums3=c(nums3,rep(i,3))
  }
  #最后一行
  nums4=c(rep(((2+row_num)*col_num)+1,col_num*3))
  nums=c(nums1,nums2,nums3,nums4)
  showLayout <- F 
  
  layout(matrix(nums,byrow = T,nrow = row_nums))
  if(showLayout) {
    layout.show(n = ((2+row_num)*col_num)+1) # 直观展示画布分布
  }
  par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2) 
  plot(1,1,
       xlab = "",xaxt = "n", # 不显示x坐标轴
       ylab = "",yaxt = "n") # 不显示y坐标轴
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white") # 背景涂黑
  text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
       (par("usr")[3]+par("usr")[4])/2,
       gname,cex = 2, col = "black") # 显示图标题
  #标题
  for (i in 1:length(vname)){
    plot(1,1,
         xlab = "",xaxt = "n", 
         ylab = "",yaxt = "n") 
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white") 
    text((par("usr")[1]+par("usr")[2])/2, 
         (par("usr")[3]+par("usr")[4])/2,
         vname[i],cex = 2, col = "black") 
  }
  #抬头和扇形图
  for (i in 1:length(g_n)){
    plot(1,1,
         xlab = "",xaxt = "n", # 不显示x坐标轴
         ylab = "",yaxt = "n") # 不显示y坐标轴
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white") # 背景涂黑
    text((par("usr")[1]+par("usr")[2])/2,
         (par("usr")[3]+par("usr")[4])/2,
         paste0(g_n[i],"\n(n = ",as.numeric(table(dat[,gname])[i]),")"),cex = 2, col = "black") 
    for (j in 1:length(vname)){
      aa=as.character(unique(dat[,vname[j]]))
      pie.gname1=pie.gname[pie.gname$Var1==g_n[i],]
      pie.gname1=pie.gname[pie.gname$Var1==g_n[i] & pie.gname$Var2 %in% aa,]
      pie(pie.gname1$Pct, 
          col = vname_col[[j]], 
          border = "white",  
          radius = 1, 
          labels = NA,
          init.angle = 90)
      symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    }
    
  }
  
  plot(1,1,
       xlab = "",xaxt = "n",
       ylab = "",yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white") 
  text((par("usr")[1]+par("usr")[2])/2,
       (par("usr")[3]+par("usr")[4])/2,
       'chisq.test',cex = 2, col = "black") 
  for(i in vname){
    plot(1,1,col = "white",
         xlab = "",xaxt = "n", # 不显示x坐标轴
         ylab = "",yaxt = "n") # 不显示y坐标轴
    text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
         (par("usr")[3]+par("usr")[4])/2,
         paste0("p = ",fisher.p[i]),cex = 1.5, col = "black") # 显示图标题
    abline(h = par("usr")[3], col = "black")
  }
  plot(0,0,col = "white",
       xlab = "",xaxt = "n", # 不显示x坐标轴
       ylab = "",yaxt = "n") # 不显示y坐标轴
  leg_nam=c()
  for (ii in 1:length(vname)){
    aa=as.character(unique(dat[,vname[ii]]))
    pie.gname1=pie.gname[pie.gname$Var2 %in% aa,]
    leg_nam=c(leg_nam,as.character(unique(pie.gname1$Var2)))
  }
  vname_col_name=unlist(vname_col)
  
  legend("topleft",ncol=ceiling(length(leg_nam)/2)+1,
         legend = leg_nam,
         fill = vname_col_name,
         border = NA,
         bty = "n", 
         cex = 1.2,
         x.intersp = 0.05,
         y.intersp = 1,
         text.width = 0.1, 
         horiz = F)
}


pdf('results/ana/06.subtype.cli/Fig5.pdf',height = 7,width = 15)
pie_compare_plot(dat = tcga_subtype.cli[,c("Cluster","Gender","T.Stage","N.Stage","M.Stage","Stage","Age1",'Event')],
                 gname = 'Cluster',
                 group_cols = c(ggsci::pal_jama()(9)))
dev.off()
#7、突变#####
tcga.mut.dat <- getTCGAMAFByCode('LUAD')
tcga.mut.dat <- as.data.frame(tcga.mut.dat@data)
tcga.mut.dat <- tcga.mut.dat[, c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")]

tcga.mut.dat$Variant_Classification <- 1
tcga.mut.dat <- reshape2::dcast(data = tcga.mut.dat, Hugo_Symbol ~ Tumor_Sample_Barcode)
class(tcga.mut.dat)
rownames(tcga.mut.dat) <- tcga.mut.dat$Hugo_Symbol
tcga.mut.dat <- tcga.mut.dat[, -1]

colnames(tcga.mut.dat) <- paste0(colnames(tcga.mut.dat), '-01')
mut.samples <- intersect(colnames(tcga.mut.dat), tcga_subtype.cli$Samples)


tcga.mut.dat <- tcga.mut.dat[, mut.samples]
tcga_mut_cli <- tcga_subtype.cli[mut.samples, ]

tcga.mut.dat.freq <- as.data.frame(rowSums(tcga.mut.dat))
colnames(tcga.mut.dat.freq) <- 'Freq'
tcga.mut.dat.freq$Genes <- rownames(tcga.mut.dat.freq)
library(dplyr)
str(tcga.mut.dat.freq)
head(tcga.mut.dat.freq)
tcga.mut.dat.freq <- dplyr::arrange(tcga.mut.dat.freq, desc(Freq))
head(tcga.mut.dat.freq)
dim(tcga.mut.dat.freq)
write.csv(tcga.mut.dat.freq,
          file = 'results/ana/06.subtype.cli/tcga.subtype.mut.gene.csv')

write.csv(tcga.mut.dat.freq,
          file = 'results/files/tcga.subtype.mut.gene.csv')

mut.genes <- rownames(tcga.mut.dat.freq)[tcga.mut.dat.freq$Freq > 3]
length(mut.genes)

tcga.mut.dat <- ifelse(tcga.mut.dat > 0, 'Mutant', 'WildType')

dim(tcga.mut.dat)


mut.res <- data.frame(clust1 = NA,
                      clust2 = NA,
                      clust3=NA)
mut.p <- c()
for (ge in mut.genes) {
  print(ge)
  tmp <- table(tcga.mut.dat[ge, ], tcga_mut_cli$Cluster)
  pvalue <- fisher.test(tmp)
  mut.p <- c(mut.p, pvalue$p.value)
  mut.res <- rbind(mut.res, tmp[1, ])
}
mut.res <- na.omit(mut.res)
rownames(mut.res) <- mut.genes
class(mut.res)
mut.res$P.value <- mut.p

table(mut.res$P.value < 0.05)
mut.res.filtered <- mut.res[which(mut.res$P.value < 0.05), ]
mut.res.filtered
dim(mut.res.filtered)
#847   4
write.csv(mut.res.filtered,
          file = 'results/ana/06.subtype.cli/tcga.subtype.filter.genes.csv')

write.csv(mut.res.filtered,
          file = 'results/files/tcga.subtype.filter.genes.csv')

mut.plot.dat <- tcga.mut.dat[rownames(mut.res.filtered), ]
mut.plot.dat[mut.plot.dat == 'WildType'] <- NA
# rownames(mut.plot.dat) <- paste(rownames(mut.plot.dat), signif(mut.res.filtered$P.value[1:15], 1), sep = ' ')
library(scales)
library(ggsci)
library(ComplexHeatmap)
alter_graphic <- function (graphic = c("rect", "point"), width = 1, 
                           height = 1, horiz_margin = unit(1, "pt"), 
                           vertical_margin = unit(1,  "pt"), fill = "red", col = NA, pch = 16, 
                           ...) {
  graphic = match.arg(graphic)[1]
  if (graphic == "rect") {
    if (!is.numeric(width)) {
      stop_wrap("`width` should be nummeric.")
    }
    if (!is.numeric(height)) {
      stop_wrap("`height` should be nummeric.")
    }
    if (width != 1) {
      if (missing(horiz_margin)) {
        horiz_margin = unit(0, "pt")
      }
    }
    if (height != 1) {
      if (missing(vertical_margin)) {
        vertical_margin = unit(0, "pt")
      }
    }
    fun = function(x, y, w, h) {
      w = w * width
      h = h * height
      grid.rect(x, y, w - horiz_margin * 2, h - vertical_margin * 
                  2, gp = gpar(fill = fill, col = col, ...))
    }
  }
  else if (graphic == "point") {
    fun = function(x, y, w, h) {
      grid.points(x, y, pch = pch, gp = gpar(fill = fill, 
                                             col = col, ...))
    }
  }
  return(fun)
}

# 定义突变颜色 
col = c("Mutant" = 'Black')

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  
  Mutant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Mutant"], col = NA))
  }
)

alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),	
  Mutant = alter_graphic("rect", fill = col["Mutant"])
)
heatmap_legend_param = list(title = "Alterations", 
                            at = c("Mutant"), 
                            labels = c("Mutant"))

color_cluster =  ggsci::pal_nejm()(9)[2:4]
names(color_cluster) = c('clust1', 'clust2','clust3')

c1.tcga_mut_cli <- tcga_mut_cli[tcga_mut_cli$Cluster == 'clust1', ]
c1.mut.plot.dat <- mut.plot.dat[1:20, c1.tcga_mut_cli$Samples]
# pdf('PDFs/c1.mut.plot.dat.pdf', width = 8, height = 6, onefile = F)
table(tcga_mut_cli$Cluster)
dim(c1.mut.plot.dat)
mut1 <- oncoPrint(as.matrix(c1.mut.plot.dat),
                  row_order = rownames(c1.mut.plot.dat),
                  # column_order = tcga_mut_cli$Samples,
                  alter_fun = alter_fun, 
                  col = col, 
                  # column_title = "", 
                  heatmap_legend_param = heatmap_legend_param,
                  bottom_annotation = HeatmapAnnotation(Cluster=c1.tcga_mut_cli[, c("Cluster")],
                                                        col=list("Cluster"=color_cluster),
                                                        show_annotation_name = TRUE,
                                                        gap=unit(1, "mm"),
                                                        na_col="grey"),
                  pct_side = "right", row_names_side = "left")
mut1
dev.off()


c2.tcga_mut_cli <- tcga_mut_cli[tcga_mut_cli$Cluster == 'clust2', ]
c2.mut.plot.dat <- mut.plot.dat[1:20, c2.tcga_mut_cli$Samples]

table(tcga_mut_cli$Cluster)
mut2 <- oncoPrint(as.matrix(c2.mut.plot.dat),
                  row_order = rownames(c2.mut.plot.dat),
                  # column_order = tcga_mut_cli$Samples,
                  alter_fun = alter_fun, 
                  col = col, 
                  # column_title = "", 
                  heatmap_legend_param = heatmap_legend_param,
                  bottom_annotation = HeatmapAnnotation(Cluster=c2.tcga_mut_cli[, c("Cluster")],
                                                        col=list("Cluster"=color_cluster),
                                                        show_annotation_name = TRUE,
                                                        gap=unit(1, "mm"),
                                                        na_col="grey"),
                  pct_side = "right", row_names_side = "left")
mut2
dev.off()

c3.tcga_mut_cli <- tcga_mut_cli[tcga_mut_cli$Cluster == 'clust3', ]
c3.mut.plot.dat <- mut.plot.dat[1:20, c3.tcga_mut_cli$Samples]

table(tcga_mut_cli$Cluster)
mut3 <- oncoPrint(as.matrix(c3.mut.plot.dat),
                  row_order = rownames(c3.mut.plot.dat),
                  # column_order = tcga_mut_cli$Samples,
                  alter_fun = alter_fun, 
                  col = col, 
                  # column_title = "", 
                  heatmap_legend_param = heatmap_legend_param,
                  bottom_annotation = HeatmapAnnotation(Cluster=c3.tcga_mut_cli[, c("Cluster")],
                                                        col=list("Cluster"=color_cluster),
                                                        show_annotation_name = TRUE,
                                                        gap=unit(1, "mm"),
                                                        na_col="grey"),
                  pct_side = "right", row_names_side = "left")
mut3
dev.off()
pdf('results/ana/06.subtype.cli/Fig5B.pdf', width = 12, height = 5, onefile = F)
mut1 + mut2 +mut3
dev.off()


tcga.immu.lands.p1=readMatrix(paste0(MG_Grobal_baseFolder,'/source/PMC5982584_supplement_2.txt'))
colnames(tcga.immu.lands.p1)
tcga.immu.lands.p1<-tcga.immu.lands.p1[tcga.immu.lands.p1$`TCGA Study`=='LUAD',]
rownames(tcga.immu.lands.p1)=paste0(rownames(tcga.immu.lands.p1),'-01')

tcga.subtype_forAlt=tcga.subtype
rownames(tcga.subtype_forAlt)=substr(rownames(tcga.subtype_forAlt),1,15)

table(is.na(match(row.names(tcga.subtype_forAlt),row.names(tcga.immu.lands.p1))))

tcga.immu.lands.p1=tcga.immu.lands.p1[intersect(row.names(tcga.subtype_forAlt),row.names(tcga.immu.lands.p1)),]
dim(tcga.immu.lands.p1)
colnames(tcga.immu.lands.p1)

table(tcga.immu.lands.p1$`TCGA Subtype`,tcga.subtype_forAlt[rownames(tcga.immu.lands.p1),1])
plotMutiBar(table(tcga.immu.lands.p1$`TCGA Subtype`
                  ,tcga.subtype_forAlt[rownames(tcga.immu.lands.p1),1]))

table(tcga.immu.lands.p1$`Immune Subtype`,tcga.subtype_forAlt[rownames(tcga.immu.lands.p1),1])
plotMutiBar(table(tcga.immu.lands.p1$`Immune Subtype`
                  ,tcga.subtype_forAlt[rownames(tcga.immu.lands.p1),1]))

colnames(tcga.immu.lands.p1)

col.selected=c('Homologous Recombination Defects','Fraction Altered','Number of Segments','Nonsilent Mutation Rate')
length(col.selected)


fig4b=list()
fig4b[[1]]=mg_violin(data.frame(tcga.subtype_forAlt[rownames(tcga.immu.lands.p1),"Cluster"]
                                ,tcga.immu.lands.p1[,'Homologous Recombination Defects'])
                     ,melt = T
                     ,ylab = 'Homologous Recombination Defects'
                     # ,ylim = c(0,100)
                     ,test_method = 'wilcox.test'
                     ,cmp_test_method = 'wilcox.test'
                     ,legend.pos = 'tr'
                     ,jitter=T
                     ,show_compare = T)
fig4b[[1]]
fig4b[[2]]=mg_violin(data.frame(tcga.subtype_forAlt[rownames(tcga.immu.lands.p1),"Cluster"]
                                ,tcga.immu.lands.p1[,'Fraction Altered'])
                     ,melt = T
                     ,ylab = 'Fraction Altered'
                     # ,ylim = c(0,100)
                     ,test_method = 'wilcox.test'
                     ,cmp_test_method = 'wilcox.test'
                     ,legend.pos = 'tr'
                     ,jitter=T
                     ,show_compare = T)
fig4b[[2]]
fig4b[[3]]=mg_violin(data.frame(tcga.subtype_forAlt[rownames(tcga.immu.lands.p1),"Cluster"]
                                ,tcga.immu.lands.p1[,'Number of Segments'])
                     ,melt = T
                     ,ylab = 'Number of Segments'
                     # ,ylim = c(0,100)
                     ,test_method = 'wilcox.test'
                     ,cmp_test_method = 'wilcox.test'
                     ,legend.pos = 'tr'
                     ,jitter=T
                     ,show_compare = T)
fig4b[[3]]

fig4b[[4]]=mg_violin(na.omit(data.frame(tcga.subtype_forAlt[rownames(tcga.immu.lands.p1),"Cluster"]
                                        ,tcga.immu.lands.p1[,'Nonsilent Mutation Rate']))
                     ,melt = T
                     ,ylab = 'Tumor mutation burden'
                     #,ylim = c(0,40)
                     ,test_method = 'wilcox.test'
                     ,cmp_test_method = 'wilcox.test'
                     ,legend.pos = 'tr'
                     ,jitter=T
                     ,show_compare = T)
fig4b[[4]]
fig4b<-mg_merge_plot(fig4b,nrow = 1,ncol = 4,common.legend = T)

ggsave('results/ana/06.subtype.cli/Fig5c.pdf',fig4b,height = 5,width = 12)
#8、GSEA 功能富集分析#####
dir.create('results/ana/07.GSEA')

tcga.subtype.1=tcga.subtype[,c("Samples","Cluster")]
tcga.subtype.1$Cluster1=ifelse(tcga.subtype.1$Cluster=='clust1','clust1','no_clust1')
tcga.subtype.1$Cluster2=ifelse(tcga.subtype.1$Cluster=='clust2','clust2','no_clust2')
tcga.subtype.1$Cluster3=ifelse(tcga.subtype.1$Cluster=='clust3','clust3','no_clust3')
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  table(degs_C1_C3$DEG$adj.P.Val<0.05)
  table(degs_C1_C3$DEG$adj.P.Val<0.01)
  
  ## 未过滤任何基因
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  # degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  # degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  # 
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  # gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$SYMBOL 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 降序排列
  return(geneList)
}
dim(tcga_tmp_log2_T)
table(tcga.subtype.1$Cluster)
tcga.geneList1=getGeneFC(gene.exp=tcga_tmp_log2_T[,tcga.subtype.1$Samples],
                        group=tcga.subtype.1$Cluster1 ,
                        ulab='clust1',
                        dlab='no_clust1')
tcga.geneList2=getGeneFC(gene.exp=tcga_tmp_log2_T[,tcga.subtype.1$Samples],
                         group=tcga.subtype.1$Cluster2 ,
                         ulab='clust2',
                         dlab='no_clust2')
tcga.geneList3=getGeneFC(gene.exp=tcga_tmp_log2_T[,tcga.subtype.1$Samples],
                         group=tcga.subtype.1$Cluster3 ,
                         ulab='clust3',
                         dlab='no_clust3')

h.all.gmt<-read.gmt("/pub1/data/mg_projects/users/wangtl/public/GSEA.gmt/h.all.v7.5.1.symbols.gmt")

tcga.hallmark.gsea1<-GSEA(tcga.geneList1,TERM2GENE = h.all.gmt,seed=T)
tcga.hallmark.gsea2<-GSEA(tcga.geneList2,TERM2GENE = h.all.gmt,seed=T)
tcga.hallmark.gsea3<-GSEA(tcga.geneList3,TERM2GENE = h.all.gmt,seed=T)

########## 热图展示 异常的通路
tcga.hallmark.gsea1.res=tcga.hallmark.gsea1@result
tcga.hallmark.gsea2.res=tcga.hallmark.gsea2@result
tcga.hallmark.gsea3.res=tcga.hallmark.gsea3@result


table(tcga.hallmark.gsea1.res$NES>0)
# FALSE 
# 18
table(tcga.hallmark.gsea2.res$NES>0)
# TRUE 
# 23 
table(tcga.hallmark.gsea3.res$NES>0)
# TRUE 
# 16 
rownames(tcga.hallmark.gsea1.res)=gsub("HALLMARK_","",rownames(tcga.hallmark.gsea1.res))
rownames(tcga.hallmark.gsea2.res)=gsub("HALLMARK_","",rownames(tcga.hallmark.gsea2.res))
rownames(tcga.hallmark.gsea3.res)=gsub("HALLMARK_","",rownames(tcga.hallmark.gsea3.res))

write.table(tcga.hallmark.gsea1.res,'results/ana/07.GSEA/tcga.hallmark.gsea.res.c1vsno_c1.txt',quote = F,row.names = F,sep='\t')
write.table(tcga.hallmark.gsea2.res,'results/ana/07.GSEA/tcga.hallmark.gsea.res.c2vsno_c2.txt',quote = F,row.names = F,sep='\t')
write.table(tcga.hallmark.gsea3.res,'results/ana/07.GSEA/tcga.hallmark.gsea.res.c3vsno_c3.txt',quote = F,row.names = F,sep='\t')



hallmark.union=Reduce(union,list(rownames(tcga.hallmark.gsea1.res),
                                 rownames(tcga.hallmark.gsea2.res),
                                 rownames(tcga.hallmark.gsea3.res)))
length(hallmark.union)

hallmark.heatmap.dat=matrix(0,nrow = 3,ncol = length(hallmark.union))
rownames(hallmark.heatmap.dat)=c('clust1 vs no_clust1', 'clust2 vs no_clust2','clust3 vs no_clust3')
colnames(hallmark.heatmap.dat)=hallmark.union

hallmark.heatmap.dat[1,match(rownames(tcga.hallmark.gsea1.res),colnames(hallmark.heatmap.dat))]=tcga.hallmark.gsea1.res$NES

hallmark.heatmap.dat[2,match(rownames(tcga.hallmark.gsea2.res),colnames(hallmark.heatmap.dat))]=tcga.hallmark.gsea2.res$NES
hallmark.heatmap.dat[3,match(rownames(tcga.hallmark.gsea3.res),colnames(hallmark.heatmap.dat))]=tcga.hallmark.gsea3.res$NES

range(hallmark.heatmap.dat)
library(ComplexHeatmap)
pdf('results/ana/07.GSEA/Fig6.pdf',height = 5,width =12)
Heatmap(as.matrix(hallmark.heatmap.dat)
        , name = "NES"
        , col = circlize::colorRamp2(c(-3, 0, 3), c('blue', 'white', 'yellow'))
        , border = TRUE
        , show_column_names = T
        , show_column_dend = F
        , show_row_dend = F
        , cluster_columns=T
        , cluster_rows=F
        , rect_gp = gpar(col = "white", lwd = 1)
        , row_names_gp = gpar(fontsize = 10)
)
dev.off()


#tcga.onc<-immu_oncogenic_pathways(exp_data = tcga_tmp_log2_T)
#save(tcga.onc,file = 'tcga.onc.RData')
load('tcga.onc.RData')
fig6d<-mg_PlotMutiBoxplot(data = tcga.onc[tcga_subtype.cli$Samples,],
                          group = tcga_subtype.cli$Cluster,
                          group_cols = ggsci::pal_nejm()(9)[2:4],
                          test_method = 'kruskal.test',
                          xangle = 30,
                          add = 'boxplot', 
                          legend.pos = 'top',
                          ylab = 'pathway ssGSEA score')
fig6d

ggsave('results/ana/07.GSEA/Fig6d.pdf',fig6d,height = 5,width = 12)
 



#9、免疫特征#####
dir.create('results/ana/08.subtype.immnu')
#tcga.cib<-immu_CIBERSORT(exp_data = tcga_tmp_log2_T)
#save(tcga.cib,file = 'tcga.cib.RData')
load('tcga.cib.RData')
fig7a<-mg_PlotMutiBoxplot(data = tcga.cib[tcga_subtype.cli$Samples,1:22],
                          group = tcga_subtype.cli$Cluster,
                          group_cols = ggsci::pal_nejm()(9)[2:4],
                          test_method = 'kruskal.test',
                          xangle = 30,
                          add = 'boxplot', 
                          legend.pos = 'top',
                          ylab = 'Score')
fig7a
#tcga.ssgsea<-immu_ssgsea(exp =tcga_tmp_log2_T,isTCGA = T )
#save(tcga.ssgsea,file='tcga.ssgsea.RData')
# load('tcga.ssgsea.RData')
# fig7a<-mg_PlotMutiBoxplot(data = tcga.ssgsea[tcga_subtype.cli$Samples,],
#                           group = tcga_subtype.cli$Cluster,
#                           group_cols = ggsci::pal_nejm()(9)[2:4],
#                           test_method = 'kruskal.test',
#                           xangle = 30,
#                           add = 'boxplot', 
#                           legend.pos = 'top',
#                           ylab = 'Score')
# fig7a
#write.table(tcga.ssgsea,'results/ana/08.subtype.immnu/tcga.ssgsea.txt',quote = F,row.names = T,sep='\t')
#tcga.esti<-immu_estimate(exp = tcga_tmp_log2_T,platform = 'illumina',isTCGA = T)
#save(tcga.esti,file = 'tcga.esti.RData')
load('tcga.esti.RData')
fig7b<-mg_PlotMutiBoxplot(data = tcga.esti[tcga_subtype.cli$Samples,1:3],
                          group = tcga_subtype.cli$Cluster,
                          group_cols = ggsci::pal_nejm()(9)[2:4],
                          test_method = 'kruskal.test',
                          xangle = 30,
                          add = 'boxplot', 
                          legend.pos = 'top',
                          ylab = 'Score')
fig7b
#免疫检查点
tcga.icg<-immu_ICGs(exp = tcga_tmp_log2_T)
fig7c<-mg_PlotMutiBoxplot(data = tcga.icg[tcga_subtype.cli$Samples,],
                          group = tcga_subtype.cli$Cluster,
                          group_cols = ggsci::pal_nejm()(9)[2:4],
                          test_method = 'kruskal.test',
                          xangle = 30,
                          add = 'boxplot', 
                          legend.pos = 'top',
                          ylab = 'Score')
fig7c
fig7ab<-mg_merge_plot(fig7a,fig7b,nrow = 1,ncol = 2,labels = c('A','B'),widths = c(3,1),align = 'h')
fig7abc<-mg_merge_plot(fig7ab,fig7c,nrow = 2,ncol = 1,labels = c('','C'))
ggsave('results/ana/08.subtype.immnu/Fig7.pdf',fig7abc,height = 12,width = 15)
#免疫相关通路
immnu_patwhays<-list.files('origin_datas/immnue_pathway/',pattern = '.gmt$')
immnu_patwhay.gene= as.data.frame(data.table::rbindlist(lapply(immnu_patwhays, function(x){
  print(x)
  a=clusterProfiler::read.gmt(paste0('origin_datas/immnue_pathway/',x))
  return(a)
})))
class(immnu_patwhay.gene)
head(immnu_patwhay.gene)
table(immnu_patwhay.gene$ont)
immnu_patwhay.gene$ont=gsub('KEGG_','',immnu_patwhay.gene$ont)
immnu_patwhay.gene=immnu_patwhay.gene[which(immnu_patwhay.gene$gene %in% rownames(tcga_tmp_log2)),]
write.table(immnu_patwhay.gene,'results/ana/08.subtype.immnu/immnu_patwhay.gene.txt',quote = F,row.names = F,sep='\t')
#差异分析筛选差异的免疫基因
immnu_patwhay.gene.exp<-tcga_tmp_log2_T[intersect(unique(immnu_patwhay.gene$gene),
                                                  rownames(tcga_tmp_log2_T)),]
diff_exp<-function(dat=immnu_patwhay.gene.exp,group=tcga.subtype$Cluster){
  gr<-as.character(unique(group))
  gg=data.frame(group=group)
  rownames(gg)=colnames(dat)
  gg1=rownames(gg[which(gg$group==gr[1]),,drop=F])
  gg2=rownames(gg[which(gg$group==gr[2]),,drop=F])
  gg3=rownames(gg[which(gg$group==gr[3]),,drop=F])
  
  p.va=data.frame()
  for (i in 1:nrow(dat)){
    dd1=as.numeric(dat[i,gg1])
    dd2=as.numeric(dat[i,gg2])
    dd3=as.numeric(dat[i,gg3])
    p.va1<-kruskal.test(list(dd1,dd2,dd3))$p.value
    if (nrow(p.va)>0){
      p.va=rbind.data.frame(p.va,
                            data.frame(gene=rownames(dat)[i],p=p.va1))
    }else{
      p.va=data.frame(gene=rownames(dat)[i],p=p.va1)
    }
    
  }
  return(p.va)
}

immnu_gene.p<-diff_exp(dat = immnu_patwhay.gene.exp[,tcga.subtype$Samples],
                       group = tcga.subtype$Cluster)
head(immnu_gene.p)
immnu_gene.p.fit<-immnu_gene.p[which(immnu_gene.p$p<0.001),]

immnu.gene<-intersect(immnu_gene.p.fit$gene,unique(immnu_patwhay.gene$gene))
length(immnu.gene)

#绘制热图
library(ComplexHeatmap)
library(circlize)

immnu_patwhay.gene.exp<-tcga_tmp_log2[immnu.gene,]
#Z-score
immnu_patwhay.gene.exp=t(apply(immnu_patwhay.gene.exp,1,function(x){return(mosaic::zscore(x))}))
immnu_patwhay.gene.exp[1:4,1:4]

tcga.subtype.group=tcga.subtype[order(tcga.subtype$Cluster),,drop=F]
col_anno <-HeatmapAnnotation(Cluster= tcga.subtype.group$Cluster)
names(col_anno)=c('Cluster')
immnu_patwhay.gene.exp=immnu_patwhay.gene.exp[,tcga.subtype.group$Samples]
#过滤
immnu_patwhay.gene=immnu_patwhay.gene[which(immnu_patwhay.gene$gene %in% immnu.gene),]
immnu_patwhay.gene1=as.data.frame(table(immnu_patwhay.gene$ont))
immnu_patwhay.gene=immnu_patwhay.gene[which(immnu_patwhay.gene$ont %in% immnu_patwhay.gene1[immnu_patwhay.gene1$Freq>1,'Var1']),]
immnu_patwhays=as.character(unique(immnu_patwhay.gene$ont))
fig7d=list()
for(i in 1:length(immnu_patwhays)){
  print(i)
  gene=immnu_patwhay.gene[which(immnu_patwhay.gene$ont==immnu_patwhays[i]),"gene"]
  mat1=as.matrix(immnu_patwhay.gene.exp[gene,])
  
  gene_ann = rowAnnotation(df = data.frame(pathway=rep(immnu_patwhays[i],length(gene)),
                                           row.names = gene))
  names(gene_ann)='pathway'
  if (i==1){
    fig7d[[i]]= Heatmap(matrix = mat1,
                        name="Row-Zscore Expression",
                        cluster_columns = F,
                        cluster_rows = T,
                        left_annotation= gene_ann,
                        col = colorRamp2(c(-3,0,3), c("blue", "white", "red")),
                        show_column_names = F,#row_title = immnu_patwhays[i],
                        show_row_names = T,show_heatmap_legend = T,
                        top_annotation = col_anno, 
                        show_row_dend = F)
  }else{
    fig7d[[i]]= Heatmap(matrix = mat1,
                        name="Row-Zscore Expression",
                        cluster_columns = F,
                        cluster_rows = T,
                        left_annotation= gene_ann,
                        col = colorRamp2(c(-3,0,3), c("blue", "white", "red")),
                        show_column_names = F,#row_title = immnu_patwhays[i],
                        show_row_names = T,show_heatmap_legend = F,
                        # top_annotation = col_anno, 
                        show_row_dend = F)
  }
  
}
length(fig7d)
# i=15
# immnu_patwhays[i]
# paste0(immnu_patwhay.gene[which(immnu_patwhay.gene$ont==immnu_patwhays[i]),"gene"],collapse = ',')

pdf('results/ana/08.subtype.immnu/Fig7D.pdf',height = 20,width = 15)
draw(fig7d[[1]] %v% fig7d[[2]] %v% fig7d[[3]] %v% fig7d[[4]] %v% fig7d[[5]] %v% fig7d[[6]] %v% fig7d[[7]] %v% fig7d[[8]] %v% fig7d[[9]] %v% fig7d[[10]] %v% fig7d[[11]] %v% fig7d[[12]] %v% fig7d[[13]] %v% fig7d[[14]]%v% fig7d[[15]],
     heatmap_legend_side = "bottom", # 热图注释放底部
     annotation_legend_side = "bottom") # 顶部注释放底部) 
dev.off()
#10、风险模型的构建#####
dir.create('results/ana/09.module')
tcga.diff1<-mg_limma_DEG(exp = tcga_tmp_log2_T[,tcga.subtype.1$Samples],
                         group = tcga.subtype.1$Cluster1,
                         ulab = 'clust1',dlab = 'no_clust1')
tcga.diff2<-mg_limma_DEG(exp = tcga_tmp_log2_T[,tcga.subtype.1$Samples],
                         group = tcga.subtype.1$Cluster2,
                         ulab = 'clust2',dlab = 'no_clust2')
tcga.diff3<-mg_limma_DEG(exp = tcga_tmp_log2_T[,tcga.subtype.1$Samples],
                         group = tcga.subtype.1$Cluster3,
                         ulab = 'clust3',dlab = 'no_clust3')
fc_fit=log2(1.5)
p_fit=0.05

tcga.diff1.fit<-tcga.diff1$DEG[which(abs(tcga.diff1$DEG$logFC)>fc_fit & tcga.diff1$DEG$adj.P.Val<p_fit),]
tcga.diff2.fit<-tcga.diff2$DEG[which(abs(tcga.diff2$DEG$logFC)>fc_fit & tcga.diff2$DEG$adj.P.Val<p_fit),]
tcga.diff3.fit<-tcga.diff3$DEG[which(abs(tcga.diff3$DEG$logFC)>fc_fit & tcga.diff3$DEG$adj.P.Val<p_fit),]

write.table(tcga.diff1$DEG,'results/ana/09.module/tcga.diff.clust1.txt',quote = F,row.names = T,sep='\t')
write.table(tcga.diff1$DEG,'results/files/tcga.diff.clust1.txt',quote = F,row.names = T,sep='\t')

write.table(tcga.diff2$DEG,'results/ana/09.module/tcga.diff.clust2.txt',quote = F,row.names = T,sep='\t')
write.table(tcga.diff2$DEG,'results/files/tcga.diff.clust2.txt',quote = F,row.names = T,sep='\t')

write.table(tcga.diff3$DEG,'results/ana/09.module/tcga.diff.clust3.txt',quote = F,row.names = T,sep='\t')
write.table(tcga.diff3$DEG,'results/files/tcga.diff.clust3.txt',quote = F,row.names = T,sep='\t')
table(tcga.diff1.fit$logFC>0)
# FALSE  TRUE 
#1694   139 
table(tcga.diff2.fit$logFC>0)
# FALSE  TRUE 
#0   18
table(tcga.diff3.fit$logFC>0)
# FALSE  TRUE 
# 77  1187

all.diff.gene=unique(c(rownames(tcga.diff1.fit),
                       rownames(tcga.diff2.fit),
                       rownames(tcga.diff3.fit)))
length(all.diff.gene)
#1915
write.table(all.diff.gene,'results/ana/09.module/all.diff.gene.txt',quote = F,row.names = F,sep='\t',col.names = F)
write.table(all.diff.gene,'results/files/all.diff.gene.txt',quote = F,row.names = F,sep='\t',col.names = F)

coxRun<-function(dat){
  library(survival)
  colnames(dat)=c('time','status','AS')  
  dat=dat[which(!is.na(dat[,1])&!is.na(dat[,3])&!is.na(dat[,2])),]
  #print(nrow(dat))
  if(nrow(dat)<10){
    print(paste0('Sample Num is small:',nrow(dat)))
    return(c(NA,NA,NA,NA))
  }
  #if(quantile(dat[,3])['25%']==quantile(dat[,3])['50%']) return(c(NA,NA,NA,NA))
  fmla <- as.formula("Surv(time, status) ~AS")
  if(table(dat[,2])[1]>1&table(dat[,2])[2]>1){
    cox <- survival::coxph(fmla, data = dat)
    re=c(summary(cox)[[7]][5],summary(cox)[[7]][2],summary(cox)[[8]][3],summary(cox)[[8]][4])
    return(re)
  }else{
    return(c(NA,NA,NA,NA))
  }
}
cox_batch<-function(dat,time,event){
  t.inds=which(!is.na(time)&!is.na(event))
  dat1=dat[,t.inds]
  os=time[t.inds]
  ev=event[t.inds]
  
  ct=sum(ev%in%c(0,1))
  if(ct!=length(ev)){
    print('event must be 0(alive) or 1(dead)')
    return(NULL)
  }
  
  res=t(apply(dat1, 1, function(x){
    ep=as.numeric(as.character(x))
    ind2=which(!is.na(ep))
    print(length(ind2))
    if(length(ind2)>1){
      os1=os[ind2]
      ev1=ev[ind2]
      ep1=ep[ind2]
      return(coxRun(data.frame(os1,ev1,ep1)))
    }else{
      return(c(NA,NA,NA,NA))
    }
  }))
  colnames(res)=c('p.value','HR','Low 95%CI','High 95%CI')
  row.names(res)=row.names(dat1)
  return(as.data.frame(res))
}

mg_risksocre.sig<-function(dat,os,os.time){
  #dat=tcga_cellage.exp_for
  #os=tcga_cli$OS
  #os.time=tcga_cli$OS.time
  tcga_dat_m<-cbind.data.frame(os=os,os.time=os.time,dat)
  tcga.cox <- cox_batch(t(tcga_dat_m[,3:ncol(tcga_dat_m)]),
                        time = tcga_dat_m$os.time,
                        event = tcga_dat_m$os)
  
  return(tcga.cox)
}

get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  options(ggrepel.max.overlaps=Inf)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10,
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=sig.coef,lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=cox1$coefficients,model=mult_results))
}

#单因素cox分析
tcga_tmp_log2.1=tcga_tmp_log2_T
rownames(tcga_tmp_log2.1)=gsub('-','__',rownames(tcga_tmp_log2.1))
all.diff.gene.1=gsub('-','__',all.diff.gene)

tcga_diff.exp<-t(tcga_tmp_log2.1[all.diff.gene.1,tcga_subtype.cli$Samples])
tcga.cox <- mg_risksocre.sig(dat = tcga_diff.exp,
                             os = tcga_subtype.cli$OS,
                             os.time = tcga_subtype.cli$OS.time)


cox.p.fit=0.0001
tcga.cox.fit<-tcga.cox[tcga.cox$p.value<cox.p.fit,]
dim(tcga.cox.fit)#23
tcga.cox$coef=log(tcga.cox$HR)
tcga.cox$Gene=rownames(tcga.cox)
tcga.cox$type=rep('None',nrow(tcga.cox))
tcga.cox$type[which(tcga.cox$p.value<cox.p.fit & tcga.cox$coef>0)]='Risk'
tcga.cox$type[which(tcga.cox$p.value<cox.p.fit & tcga.cox$coef<0)]='Protective'
fig8a <- ggplot(data = tcga.cox,
                aes(x = coef,
                    y = -log10(p.value)))+
  geom_point(alpha=0.4, size=3.5, aes(color=type))+
  scale_color_manual(values=c(mg_colors[2],'grey',mg_colors[1]),
                     limits = c("Protective",'None', "Risk"),name='State')+
  geom_hline(yintercept = -log10(cox.p.fit),lty=4,col="black",lwd=0.8)+
  ylab('-log10(pvalue)')+xlab('Cox coefficient')+
  theme_bw()+
  theme(axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
        legend.text=element_text(face="plain", family="Times", colour="black"),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black" ),#设置图例的总标题的字体属性
        legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill = NA, colour = NA)
  )
# +ggrepel::geom_text_repel(data=module.genes.cox_use[which(module.genes.cox_use$p.value<0.05),],aes(label=Gene))
fig8a
table(tcga.cox$type)
# None Protective       Risk 
#1892         14          9
write.table(tcga.cox,'results/ana/09.module/tcga.cox.txt',quote = F,row.names = F,sep='\t')
write.table(tcga.cox,'results/files/tcga.cox.txt',quote = F,row.names = F,sep='\t')
tcga.lasso<-get_riskscore.lasso(dat = tcga_diff.exp[,rownames(tcga.cox.fit)],
                                os=tcga_subtype.cli$OS,
                                os.time = tcga_subtype.cli$OS.time,
                                labels=c('B','C'))

length(names(tcga.lasso$lasso.gene))#11
tcga.lasso$lambda.min
#0.02003762
tcga.module.risk<-get_riskscore(dat=tcga_diff.exp[,names(tcga.lasso$lasso.gene)],
  #dat=tcga_diff.exp[,rownames(tcga.cox.fit)],
  os=tcga_subtype.cli$OS,
  os.time=tcga_subtype.cli$OS.time,
  step=F,
  direction=c("both", "backward", "forward")[1])

length(names(tcga.module.risk$module.gene))#11
names(tcga.module.risk$module.gene)
# "CD200R1" "CLEC17A" "ZC3H12D" "GNG7"    "SNX30"   "CDCP1"   "NEIL3"   "IGF2BP1" "RHOV"    "ABCC2"   "KRT81" 
tcga.module.risk$model
#-0.391*CD200R1+-0.141*CLEC17A+-0.076*ZC3H12D+-0.049*GNG7+-0.038*SNX30+0.229*CDCP1+0.05*NEIL3+0.05*IGF2BP1+0.081*RHOV+0.006*ABCC2+0.073*KRT81"
plotRiskScoreModel_1=function(riskScore,dat,time,event,cutoff,hetTitle='z-score of expression',hetColor=c('green','black','red')){
  srt.inds=order(riskScore)
  dat=dat[srt.inds,]
  time=time[srt.inds]
  event=event[srt.inds]
  riskScore=riskScore[srt.inds]
  library(ggplot2)
  dt1=data.frame(V1=1:length(riskScore),V2=riskScore,RiskType=ifelse(riskScore>cutoff,'High','Low')) 
  p1=ggplot(dt1, aes(x = V1, y = V2, colour = RiskType,fill=RiskType)) +
    geom_bar(stat = 'identity', position = 'dodge')+
    ggsci::scale_fill_simpsons()+ggsci::scale_color_simpsons()+theme_bw()
  p1=p1+ylab('RiskScore')+
    theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_blank()
          ,axis.title.x=element_blank(),legend.position=c(1,0), legend.justification=c(1,0)
          ,legend.background = element_rect(fill = NA, colour = NA)
          ,plot.margin=unit(c(0.1, 0.1, 0, 0.1), "inches")
          ,legend.title = element_text(family="Times",face="plain")
          ,legend.text = element_text(family="Times",face="plain"))
  dt2=data.frame(V1=1:length(riskScore),V2=time,Status=ifelse(event==1,'Dead','Alive'))  
  p2=ggplot(dt2, aes(x = V1, y = V2, colour = Status,shape =Status)) +
    geom_point()+ggsci::scale_color_simpsons()+theme_bw()
  p2=p2+ylab('Time')+
    theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_blank()
          ,axis.title.x=element_blank(),legend.position=c(1,1), legend.justification=c(1,1)
          ,legend.background = element_rect(fill = NA, colour = NA)
          ,plot.margin=unit(c(0, 0.1, 0, 0.1), "inches")
          ,legend.title = element_text(family="Times",face="plain")
          ,legend.text = element_text(family="Times",face="plain"))
  
  
  data=as.data.frame(scale(dat))
  hc.r = hclust(dist(t(data)))
  data=data[,hc.r$order]
  data$ID <- 1:nrow(dat)
  #colnames(data)
  data_m <- reshape2::melt(data, id.vars=c("ID"))
  colnames(data_m)=c('ID','V1','V2')
  data_m$V2[which(data_m$V2>mean(data_m$V2)+3*sd(data_m$V2))]=mean(data_m$V2)+3*sd(data_m$V2)
  data_m$V2[which(data_m$V2<mean(data_m$V2)-3*sd(data_m$V2))]=mean(data_m$V2)-3*sd(data_m$V2)
  
  data_m$V1=mg_str_outline(data_m$V1,isCut = T,n=50)
  #print(data_m$V1)
  #print(head(data_m))
  #data_m[1:20,]
  p3 <- ggplot(data_m, aes(x=ID,y=V1)) 
  p3 <- p3 + geom_tile(aes(fill=V2))
  p3=p3+scale_fill_gradient2(low = hetColor[1],mid=hetColor[2], high = hetColor[3])
  p3=p3+theme_bw()
  p3=p3+ labs(fill=hetTitle) 
  #p3=p3+guides(fill=guide_legend(title="New Legend Title"))
  p3=p3+xlab('Samples')
  #p3 <- p3 + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  p3=p3+theme(axis.text.y=element_text(family="Times",face="plain")
              ,axis.text.x=element_blank()
              #,axis.title.x=element_blank()
              ,axis.title.y=element_blank()
              ,legend.position='bottom'
              #,legend.justification=c(1,1)
              #,legend.background = element_rect(fill = NA, colour = NA)
              ,plot.margin=unit(c(0, 0.1, 0.1, 0.1), "inches")
              ,legend.title = element_text(family="Times",face="plain")
              ,legend.text = element_text(family="Times",face="plain"))
  
  g1=ggpubr::ggarrange(p1,p2,p3, ncol = 1, nrow = 3,heights = c(0.5,0.5,1),align = "v")
  return(g1)
}
ggplotTimeROC_1=function(time,status,score,mks=c(1,3,5)){
  #time=g.os
  #status=g.ev
  #score=as.numeric(cpm.score)
  #cx=coxRun(data.frame(time,status,score))
  #if(cx[1]<=1){
  #  score=-1*score
  #}
  roc.tm=mg_surv_pROC(time,status,score,mks)
  print('roc.tm')
  print((roc.tm))
  library(survival)
  library(ggplot2)
  mks=mg_predict_time_ymd(time,mks)
  print(mks)  
  ROC.DSST=timeROC::timeROC(T=time,
                            delta=status
                            ,marker=score,
                            cause=1,weighting="marginal",
                            times=mks,
                            iid=TRUE)
  print(ROC.DSST)
  mks=mks[which(!is.na(ROC.DSST$AUC)&ROC.DSST$AUC>0)]
  print(mks)
  if(length(mks)>0){
    if(max(ROC.DSST$AUC)<0.5){
      score=-1*score
    }
    ROC.DSST=timeROC::timeROC(T=time,
                              delta=status
                              ,marker=score,
                              cause=1,weighting="marginal",
                              times=mks,
                              iid=TRUE)
    print(ROC.DSST$times)
    if(max(ROC.DSST$times)<20){
      lb=paste0(ROC.DSST$times,'-Years')
    }else if(max(ROC.DSST$times)<365){
      lb=paste0(round(ROC.DSST$times/12,0),'-Years')
    }else{
      lb=paste0(round(ROC.DSST$times/365,0),'-Years')
    }
    
    lbs=paste0(lb,',AUC=',round(ROC.DSST$AUC,2),',95%CI(',paste0(round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,1]/100,2),'-',
                                                                 round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,2]/100,2)),')')
    #roc.tm=ROC.DSST$times[which(ROC.DSST$times>0)]
    
    #p.dat=rbind()
    #for(i in which(ROC.DSST$times>0)){
    #los=lowess(ROC.DSST$FP[,i], y=ROC.DSST$TP[,i], f = 1/3, iter = 100)
    #los$x=c(0,los$x,1)
    #los$y=c(0,los$y,1)
    # p.dat=rbind(p.dat,data.frame(los$x, y=los$y,rep(lbs[i],length(los$y)),stringsAsFactors = F))
    #}
    
    p.dat=rbind()
    print(length(roc.tm))
    for(i in 1:length(roc.tm)){
      #print(i)
      r1=roc.tm[[i]]
      x1=1-r1$specificities
      y1=r1$sensitivities
      #print(cbind(1-r1$specificities,r1$sensitivities))
      nx1=unique(x1)
      ny1=c()
      for(x in unique(x1)){
        x.inds=which(x1==x)
        if(length(x.inds)>0&x<0.5){
          ny1=c(ny1,min(y1[x.inds]))
        }else if(length(x.inds)>0){
          ny1=c(ny1,max(y1[x.inds]))
        }else{
          ny1=c(ny1,y1[x.inds][1])
        }
      }
      #print(cbind(nx1,ny1))
      p.dat=rbind(p.dat,data.frame(x=nx1, y=ny1,rep(lbs[i],length(nx1)),stringsAsFactors = F))
    }
    colnames(p.dat)=c('V1','V2','Type')
    p.dat=as.data.frame(p.dat)
    
    p1=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))
    p1=p1+geom_line(aes(colour=Type),lwd=1.1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction') +ggsci::scale_fill_simpsons()+ggsci::scale_color_simpsons()
    #p1=p1+stat_smooth(aes(colour=Type),se = FALSE, size = 1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction') 
    
    p1=p1+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_text(family="Times",face="plain")
                ,axis.title.x=element_text(family="Times",face="plain"),axis.title.y=element_text(family="Times",face="plain")
                ,plot.title=element_blank()
                ,plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches")
                ,legend.position=c(1,0)
                ,legend.justification=c(1,0)
                ,legend.background = element_rect(fill = NA, colour = NA)
                ,legend.title = element_text(family="Times",face="plain")
                ,legend.text = element_text(family="Times",face="plain"))
    return(p1)
  }else{
    return(mg_getplot_bank('No data plot by ROC!'))
  }
}
mg_predict_time_ymd=function(time,mks){
  mx=quantile(time,seq(0,1,0.1),na.rm=T)['100%']
  #mx=max(time,na.rm = T)
  if(mx<25){
  }else if(mx<365){
    mks=mks*12 
  }else{
    mks=mks*365 
  }
  mks=mks[which(mks<quantile(time,seq(0,1,0.01),na.rm=T)['100%'])]
  return(mks)
}
fig8c=plotRiskScoreModel_1(riskScore = tcga.module.risk$result$riskscorez,
                           dat = tcga_diff.exp[,names(tcga.module.risk$module.gene)],
                           time = tcga_subtype.cli$OS.time,
                           event = tcga_subtype.cli$OS,
                           cutoff = 0,
                           hetTitle = 'z-score of expression',
                           hetColor = c('yellow','white','blue'))
fig8c
fig8d<-ggplotTimeROC_1(time = tcga.module.risk$result$time/365
                       ,status = tcga.module.risk$result$status
                       ,score = tcga.module.risk$result$riskscore
                       ,mks = c(1,3,5))
fig8d
fig8e<-ggplotKMCox(dat = data.frame(tcga.module.risk$result$time/365,
                                    tcga.module.risk$result$status,
                                    #tcga.module.risk$result$Risk
                                    ifelse(tcga.module.risk$result$riskscorez>=0,'H','L')),
                   title = 'TCGA',add_text = '',
                   palette = ggsci::pal_simpsons()(9)[c(5,2)],labs = c('High','Low'))
fig8e
tcga.module.risk$result$Risk=ifelse(tcga.module.risk$result$riskscorez>=0,'High','Low')
#独立数据集验证
gse30219.module.risk<-get_riskscore(dat = t(GSE30219_exp[intersect(names(tcga.module.risk$module.gene),
                                                              rownames(GSE30219_exp)),
                                                     GSE30219_cli$Samples]),
                                os =GSE30219_cli$OS,
                                os.time = GSE30219_cli$OS.time,
                                step=F,direction=c("both", "backward", "forward")[1])

fig8f<-ggplotTimeROC_1(time = gse30219.module.risk$result$time/365
                       ,status = gse30219.module.risk$result$status
                       ,score = gse30219.module.risk$result$riskscorez
                       ,mks = c(1,3,5))
fig8f
fig8g<-ggplotKMCox(dat = data.frame(gse30219.module.risk$result$time/365,
                                    gse30219.module.risk$result$status,
                                    ifelse(gse30219.module.risk$result$riskscorez>=0,'H','L')),
                   title = 'GSE30219',add_text = '',
                   palette = ggsci::pal_simpsons()(9)[c(5,2)],labs = c('High','Low'))
fig8g
gse30219.module.risk$result$Risk=ifelse(gse30219.module.risk$result$riskscorez>=0,'High','Low')
#
gse31210.module.risk<-get_riskscore(dat = t(GSE31210_exp[intersect(names(tcga.module.risk$module.gene),
                                                                   rownames(GSE31210_exp)),
                                                         GSE31210_cli$Samples]),
                                    os =GSE31210_cli$OS,
                                    os.time = GSE31210_cli$OS.time,
                                    step=F,direction=c("both", "backward", "forward")[1])

fig8h<-ggplotTimeROC_1(time = gse31210.module.risk$result$time/365
                       ,status = gse31210.module.risk$result$status
                       ,score = gse31210.module.risk$result$riskscorez
                       ,mks = c(1,3,5))
fig8h
fig8i<-ggplotKMCox(dat = data.frame(gse31210.module.risk$result$time/365,
                                    gse31210.module.risk$result$status,
                                    ifelse(gse31210.module.risk$result$riskscorez>=0,'H','L')),
                   title = 'GSE31210',add_text = '',
                   palette = ggsci::pal_simpsons()(9)[c(5,2)],labs = c('High','Low'))
fig8i
gse31210.module.risk$result$Risk=ifelse(gse31210.module.risk$result$riskscorez>=0,'High','Low')

#GSE19188
gse19188.module.risk<-get_riskscore(dat = t(GSE19188_exp[intersect(names(tcga.module.risk$module.gene),
                                                                   rownames(GSE19188_exp)),
                                                         GSE19188_cli$Samples]),
                                    os =GSE19188_cli$OS,
                                    os.time = GSE19188_cli$OS.time,
                                    step=F,direction=c("both", "backward", "forward")[1])

fig8j<-ggplotTimeROC_1(time = gse19188.module.risk$result$time/365
                       ,status = gse19188.module.risk$result$status
                       ,score = gse19188.module.risk$result$riskscorez
                       ,mks = c(1,3,5))
fig8j
fig8k<-ggplotKMCox(dat = data.frame(gse19188.module.risk$result$time/365,
                                    gse19188.module.risk$result$status,
                                    ifelse(gse19188.module.risk$result$riskscorez>=0,'H','L')),
                   title = 'GSE19188',add_text = '',
                   palette = ggsci::pal_simpsons()(9)[c(5,2)],labs = c('High','Low'))
fig8k
gse19188.module.risk$result$Risk=ifelse(gse19188.module.risk$result$riskscorez>=0,'High','Low')
#
#GSE37745
# gse37745.module.risk<-get_riskscore(dat = t(GSE37745_exp[intersect(names(tcga.module.risk$module.gene),
#                                                                    rownames(GSE37745_exp)),
#                                                          GSE37745_cli$Samples]),
#                                     os =GSE37745_cli$OS,
#                                     os.time = GSE37745_cli$OS.time,
#                                     step=F,direction=c("both", "backward", "forward")[1])
# 
# fig8m<-ggplotTimeROC_1(time = gse37745.module.risk$result$time/365
#                        ,status = gse37745.module.risk$result$status
#                        ,score = gse37745.module.risk$result$riskscorez
#                        ,mks = c(1,3,5))
# fig8m
# fig8n<-ggplotKMCox(dat = data.frame(gse37745.module.risk$result$time/365,
#                                     gse37745.module.risk$result$status,
#                                     ifelse(gse37745.module.risk$result$riskscorez>=0,'H','L')),
#                    title = 'GSE37745',add_text = '',
#                    palette = ggsci::pal_simpsons()(9),labs = c('High','Low'))
# fig8n
# gse37745.module.risk$result$Risk=ifelse(gse37745.module.risk$result$riskscorez>=0,'High','Low')

#GSE50081
gse50081.module.risk<-get_riskscore(dat = t(GSE50081_exp[intersect(names(tcga.module.risk$module.gene),
                                                                   rownames(GSE50081_exp)),
                                                         GSE50081_cli$Samples]),
                                    os =GSE50081_cli$OS,
                                    os.time = GSE50081_cli$OS.time,
                                    step=F,direction=c("both", "backward", "forward")[1])

fig8m<-ggplotTimeROC_1(time = gse50081.module.risk$result$time/365
                       ,status = gse50081.module.risk$result$status
                       ,score = gse50081.module.risk$result$riskscorez
                       ,mks = c(1,3,5))
fig8m
fig8n<-ggplotKMCox(dat = data.frame(gse50081.module.risk$result$time/365,
                                    gse50081.module.risk$result$status,
                                    ifelse(gse50081.module.risk$result$riskscorez>=0,'H','L')),
                   title = 'GSE50081',add_text = '',
                   palette = ggsci::pal_simpsons()(9)[c(5,2)],labs = c('High','Low'))
fig8n
gse50081.module.risk$result$Risk=ifelse(gse50081.module.risk$result$riskscorez>=0,'High','Low')
#fig8abcd<-mg_merge_plot(fig8a,tcga.lasso$plot,fig8c,widths = c(1,2,1),nrow = 1,ncol =3,labels = c('A','','D'))
fig8de<-mg_merge_plot(fig8d,fig8e,nrow = 2,ncol = 1,labels = c('E','F'))
fig8abc<-mg_merge_plot(fig8a,tcga.lasso$plot,widths = c(1,2),nrow = 1,ncol =2,labels = c('A',''))
fig8fg<-mg_merge_plot(fig8f,fig8g,nrow = 2,ncol = 1,labels = LETTERS[7:8])
fig8def<-mg_merge_plot(fig8c,fig8de,fig8fg,nrow = 1,ncol = 3,labels = c('D','',''))

fig8hi<-mg_merge_plot(fig8h,fig8j,fig8m,
                      fig8i,fig8k,fig8n,
                      nrow = 2,ncol = 3,
                      labels = LETTERS[9:15])
fig8=mg_merge_plot(fig8abc,fig8def,fig8hi,nrow = 3,ncol = 1,heights = c(1,2,2))
ggsave('results/ana/09.module/Fig8.pdf',fig8,height = 22,width = 15)
#11、风险模型临床表型的差异分析####
dir.create('results/ana/10.module.cli')
head(tcga_subtype.cli)
tcga_subtype.cli.risk<-merge(tcga_subtype.cli,
                             tcga.module.risk$result[,c("Samples","riskscore","riskscorez","Risk")],
                             by='Samples')
head(tcga_subtype.cli.risk)
write.table(tcga_subtype.cli.risk,'results/ana/10.module.cli/tcga_subtype.cli.risk.txt',quote = F,row.names = F,sep='\t')

table(tcga_subtype.cli$Gender)
fig9a=list()
fig9a[[1]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("Gender","riskscore")],
                        leg = 'Gender',ylab = 'RiskScore',palette = ggsci::pal_jco()(9))
fig9a[[1]]

fig9a[[2]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("T.Stage","riskscore")],
                        leg = 'T.Stage',ylab = 'RiskScore',palette = ggsci::pal_jco()(9))
fig9a[[2]]

fig9a[[3]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("M.Stage","riskscore")],
                        leg = 'M.Stage',ylab = 'RiskScore',palette = ggsci::pal_jco()(9))
fig9a[[3]]

fig9a[[4]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("N.Stage","riskscore")],
                        leg = 'N.Stage',ylab = 'RiskScore',palette = ggsci::pal_jco()(9))
fig9a[[4]]

fig9a[[5]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("Stage","riskscore")],
                        leg = 'Stage',ylab = 'RiskScore',palette = ggsci::pal_jco()(9))
fig9a[[5]]
fig9a[[6]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("Age1","riskscore")],
                        leg = 'Age',ylab = 'RiskScore',palette = ggsci::pal_jco()(9))
fig9a[[6]]
fig9a[[7]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("Smoking","riskscore")],
                        leg = 'Smoking',ylab = 'RiskScore',palette = ggsci::pal_jco()(9))
fig9a[[7]]

fig9a[[8]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("Cluster","riskscore")],
                        leg = 'Cluster',ylab = 'RiskScore',palette = ggsci::pal_jco()(9))
fig9a[[8]]

# 
# library(ComplexHeatmap)
# pal_lancet()(9)[c(7,1)]
# tcga_subtype.cli.risk=tcga_subtype.cli.risk[order(tcga_subtype.cli.risk$riskscorez),]
# tcga_subtype.cli.risk1=tcga_subtype.cli.risk[,c("Cluster","T.Stage","N.Stage","M.Stage","Stage","Event")]
# cli.colors=list()
# for(i in 1:ncol(tcga_subtype.cli.risk1)){
#   mycolor=ggsci::pal_jco()(9)
#   var=tcga_subtype.cli.risk1[i]
#   var.clas=names(table(tcga_subtype.cli.risk1[,i]))
#   var.color=mycolor[1:length(var.clas)]
#   names(var.color)=var.clas
#   cli.colors=c(cli.colors,list(var.color))
# }
# names(cli.colors)=colnames(tcga_subtype.cli.risk1)
# 
# 
# tcga.risk.barplot=columnAnnotation(RiskScore = anno_barplot(tcga_subtype.cli.risk$riskscorez,baseline =0,bar_width=1,gp=gpar(fill=ifelse(tcga_subtype.cli.risk$riskscorez>0,'#AD002AFF','#00468BFF'))),annotation_name_side ='left',height =unit(4,'inches'))
# draw(tcga.risk.barplot)
# 
# tcga.cli.heatmap=columnAnnotation(df = tcga_subtype.cli.risk1
#                                   ,annotation_name_side='left'
#                                   ,annotation_height =unit(4,'inches')
#                                   ,col = cli.colors)
# 
# ht_list=tcga.risk.barplot %v% tcga.cli.heatmap
# ht_list
# dev.off()

plot_sankey=function(df_m){
  library(ggalluvial)
  library(ggplot2)
  library(dplyr)
  corLodes=to_lodes_form(df_m, axes = 1:ncol(df_m), id = "Cohort")
  mycol <- rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
  ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
    scale_x_discrete(expand = c(0, 0)) +  
    #用aes.flow控制线调颜色，forward说明颜色和前面一致，backward说明与后面一致
    geom_flow(width = 3/10,aes.flow = "forward") + 
    geom_stratum(alpha = .9,width = 3/10) +
    scale_fill_manual(values = mycol) +
    #size = 2代表基因名字大小
    geom_text(stat = "stratum", size = 4,color="black") +
    xlab("") + ylab("") + theme_bw() + 
    theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
    theme(panel.grid =element_blank()) + 
    theme(panel.border = element_blank()) + 
    ggtitle("") + guides(fill = FALSE)   
}
tcga_subtype.cli.risk1=tcga_subtype.cli.risk[,c("Cluster","Risk","T.Stage","N.Stage","Stage",'Gender')]
fig9b<-plot_sankey(tcga_subtype.cli.risk1)
fig9b

fig9a1<-mg_merge_plot(fig9a[1:6],nrow = 2,ncol = 3)
fig9a2<-mg_merge_plot(fig9a[7:8],nrow = 1,ncol = 2,widths = c(2,1))

fig9<-mg_merge_plot(fig9a1,fig9a2,fig9b,nrow = 3,ncol = 1,labels = c('A','','B'),heights = c(2,1,1))
ggsave('results/ana/10.module.cli/Fig9.pdf',fig9,height = 15,width = 12)
#12、单因素和多因素+列线图#####
dir.create('results/ana/11.module.sigmut')
library(forestplot)
library(survcomp)
tcga_cox_datas <- tcga_subtype.cli.risk
tcga_cox_datas <- crbind2DataFrame(tcga_cox_datas)

tcga_cox_datas$RiskType=tcga_cox_datas$Risk
table(tcga_cox_datas$RiskType)
tcga_cox_datas$RiskType <- ifelse(tcga_cox_datas$RiskType == 'High', '1High', '0Low')

table(tcga_cox_datas$T.Stage)
tcga_cox_datas$T.Stage[tcga_cox_datas$T.Stage == 'T1' | tcga_cox_datas$T.Stage == 'T2'] <- 'T1+T2'
tcga_cox_datas$T.Stage[tcga_cox_datas$T.Stage == 'T3' | tcga_cox_datas$T.Stage == 'T4'] <- 'T3+T4'

table(tcga_cox_datas$N.Stage)
tcga_cox_datas$N.Stage[tcga_cox_datas$N.Stage == 'N1' | tcga_cox_datas$N.Stage == 'N2' | tcga_cox_datas$N.Stage == 'N3'] <- 'N1+N2+N3'
table(tcga_cox_datas$M.Stage)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage == 'I' | tcga_cox_datas$Stage == 'II'] <- 'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage == 'III' | tcga_cox_datas$Stage == 'IV'] <- 'III+IV'


write.table(tcga_cox_datas,'results/ana/11.module.sigmut/tcga_cox_datas.txt',quote = F,sep = '\t',row.names = F)
# 单因素分析结果 
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

A3_T_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.Stage,
                              data=tcga_cox_datas))
A3_T_sig_cox_dat <- data.frame(Names=rownames(A3_T_sig_cox[[8]]),
                               HR = round(A3_T_sig_cox[[7]][,2],3),
                               lower.95 = round(A3_T_sig_cox[[8]][,3],3),
                               upper.95 = round(A3_T_sig_cox[[8]][,4],3),
                               p.value=round(A3_T_sig_cox[[7]][,5],3))
A3_T_sig_cox_dat

A4_N_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~N.Stage,
                              data=tcga_cox_datas))
A4_N_sig_cox_dat <- data.frame(Names=rownames(A4_N_sig_cox[[8]]),
                               HR = round(A4_N_sig_cox[[7]][,2],3),
                               lower.95 = round(A4_N_sig_cox[[8]][,3],3),
                               upper.95 = round(A4_N_sig_cox[[8]][,4],3),
                               p.value=round(A4_N_sig_cox[[7]][,5],3))
A4_N_sig_cox_dat

A5_M_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~M.Stage,
                              data=tcga_cox_datas))
A5_M_sig_cox_dat <- data.frame(Names=rownames(A5_M_sig_cox[[8]]),
                               HR = round(A5_M_sig_cox[[7]][,2],3),
                               lower.95 = round(A5_M_sig_cox[[8]][,3],3),
                               upper.95 = round(A5_M_sig_cox[[8]][,4],3),
                               p.value=round(A5_M_sig_cox[[7]][,5],3))
A5_M_sig_cox_dat


Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat


RiskType_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~riskscore,
                                  data=tcga_cox_datas))
RiskType_sig_cox_dat <- data.frame(Names=rownames(RiskType_sig_cox[[8]]),
                                   HR = round(RiskType_sig_cox[[7]][,2],3),
                                   lower.95 = round(RiskType_sig_cox[[8]][,3],3),
                                   upper.95 = round(RiskType_sig_cox[[8]][,4],3),
                                   p.value=round(RiskType_sig_cox[[7]][,5],3))
RiskType_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     #A3_T_sig_cox_dat,
                     #A4_N_sig_cox_dat,
                     #A5_M_sig_cox_dat,
                     Stage_sig_cox_dat,
                     RiskType_sig_cox_dat)

data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
data.sig$Names
rownames(data.sig) <- c("Age",
                        "Gender",
                        #"T.Stage",
                        #"N.Stage",
                        #"M.Stage",
                        "Stage",
                        "RiskScore")
data.sig$Names <- rownames(data.sig)
data.sig$p.value=ifelse(data.sig$p.value==0.000,'<0.0001',data.sig$p.value)
data.sig
pdf('results/ana/11.module.sigmut/Fig10A.pdf', width = 8, height = 5,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap = 8,lineheight = 10)
dev.off()
#多因素分析结果
#muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+Gender+T.Stage+N.Stage+M.Stage+Stage+riskscore, data=tcga_cox_datas))
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+Gender+Stage+riskscore, data=tcga_cox_datas))

muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
rownames(data.muti) <- c("Age",
                         "Gender",
                         #"T.Stage",
                         #"N.Stage",
                         #"M.Stage",
                         "Stage",
                         "RiskScore")
data.muti$Names <- rownames(data.muti)
data.muti$p.value=ifelse(data.muti$p.value==0.000,'<0.0001',data.muti$p.value)
data.muti
pdf('results/ana/11.module.sigmut/Fig10B.pdf', width = 8, height = 5,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap = 8,lineheight = 10)
dev.off()
#列线图
mg_nomogram=function(clinical_riskscore,
                     os,
                     status,
                     title='Nomogram',
                     quick=T,
                     mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3#对观测2的六个指标在列线图上进行计分展示
  #,observation=pbc[2,] #也可以不展示
  #预测3年和5年的死亡风险，此处单位是day
  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE #cox回归中需要TRUE
  #              ,showP = T #是否展示统计学差异
  #              ,droplines = F#观测2示例计分是否画线
  #,colors = mg_colors[1:3] #用前面自己定义的颜色
  #,rank="decreasing") #根据统计学差异的显著性进行变量的排序
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) #展示观测的可信区间
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}

pdf('results/ana/11.module.sigmut/fig10C.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(Stage=tcga_cox_datas$Stage,
                                RiskScore=tcga_cox_datas$riskscore),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5))
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))


#AUC
library("timeROC")
tcga_cox_auc=tcga_cox_datas
tcga_cox_auc$Stage=as.numeric(as.factor(tcga_cox_auc$Stage))
tcga_cox_auc$Riskscore=as.numeric(tcga_cox_auc$riskscore)
tcga_cox_auc$Gender=as.numeric(as.factor(tcga_cox_auc$Gender))
#tcga_cox_auc$Grade=as.numeric(as.factor(tcga_cox_auc$Grade))

# tcga_cox_auc$T.Stage=as.numeric(as.factor(tcga_cox_auc$T.Stage))
# tcga_cox_auc$Age=as.numeric(as.factor(tcga_cox_auc$Age))
# tcga_cox_auc$N.Stage=as.numeric(as.factor(tcga_cox_auc$N.Stage))
# tcga_cox_auc$M.Stage=as.numeric(as.factor(tcga_cox_auc$M.Stage))


ROC.DSST.Age=timeROC(T=tcga_cox_auc$OS.time/365,
                     delta=tcga_cox_auc$OS,
                     marker=tcga_cox_auc$Age,
                     cause=1,weighting="marginal",
                     times=c(1,2,3,4,5),
                     iid=TRUE)
ROC.DSST.Gender=timeROC(T=tcga_cox_auc$OS.time/365,
                        delta=tcga_cox_auc$OS,
                        marker=tcga_cox_auc$Gender,
                        cause=1,weighting="marginal",
                        times=c(1,2,3,4,5),
                        iid=TRUE)
ROC.DSST.Stage=timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Stage,
                       cause=1,weighting="marginal",
                       times=c(1,2,3,4,5),
                       iid=TRUE)

# ROC.DSST.Grade=timeROC(T=tcga_cox_auc$OS.time/365,
#                        delta=tcga_cox_auc$OS,
#                        marker=tcga_cox_auc$Grade,
#                        cause=1,weighting="marginal",
#                        times=c(1,2,3,4,5),
#                        iid=TRUE)

ROC.DSST.Risk=timeROC(T=tcga_cox_auc$OS.time/365,
                      delta=tcga_cox_auc$OS,
                      marker=tcga_cox_auc$Riskscore,
                      cause=1,weighting="marginal",
                      times=c(1,2,3,4,5),
                      iid=TRUE)

ROC.DSST.Nomo<-timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Riskscore,
                       other_markers=as.matrix(tcga_cox_auc[,c("Stage")]),
                       cause=1,
                       weighting="cox",
                       times=c(1,2,3,4,5),
                       iid=F)
ROC.DSST.Age$AUC
ROC.DSST.Gender$AUC
ROC.DSST.Stage$AUC
ROC.DSST.Grade$AUC
ROC.DSST.Risk$AUC
ROC.DSST.Nomo$AUC

scales::show_col(mg_colors[c(1,10:12,4,5,7:9)])

pdf('results/ana/11.module.sigmut/Fig10J.pdf',height = 5,width = 6)
plotAUCcurve(ROC.DSST.Nomo,conf.int=F,col=mg_colors[1])
plotAUCcurve(ROC.DSST.Risk,conf.int=F,col=mg_colors[2],add=TRUE)
plotAUCcurve(ROC.DSST.Stage,conf.int=F,col=mg_colors[3],add=TRUE)
plotAUCcurve(ROC.DSST.Age,conf.int=F,col=mg_colors[4],add=TRUE)
plotAUCcurve(ROC.DSST.Gender,conf.int=F,col=mg_colors[5],add=TRUE)
#plotAUCcurve(ROC.DSST.Grade,conf.int=F,col=mg_colors[6],add=TRUE)

legend("topright",c("Nomogram",'RiskScore'
                    ,'Stage'
                    ,"Age","Gender")
       ,col=mg_colors[c(1:5)],lty=1,lwd=2)

dev.off()
#13、免疫特征####
dir.create('results/ana/12.module.immnu')
load('tcga.esti.RData')
head(tcga.esti)
tcga.esti.risk<-cbind.data.frame(tcga.esti[tcga_subtype.cli.risk$Samples,],
                                 RiskScore=tcga_subtype.cli.risk$riskscorez,
                                 Risk=tcga_subtype.cli.risk$Risk)
mg_cor_point=function(x,y,method='Pearson',top_col='#D55E00',right_col='#009E73'
                      ,ylab='y expression',xlab='x expression',title=NULL
                      ,marginal.type=c("histogram", "boxplot", "density", "violin", "densigram")[1]){
  library(ggstatsplot)
  dat=data.frame(X=x,Y=y)
  tp='nonparametric'
  if(method=='Pearson'|method=='pearson'){
    tp='parametric'
  }
  g1=ggscatterstats(data = dat,
                    x = X,
                    y = Y
                    ,type = tp
                    ,xfill = top_col
                    ,yfill = right_col
                    ,xlab = xlab
                    ,ylab=ylab
                    ,marginal.type = marginal.type
                    ,title = title)
  return(g1)
}

fig11a<-mg_cor_point(x = tcga.esti.risk$ImmuneScore,y = tcga.esti.risk$RiskScore,
                     xlab = 'ImmuneScore',ylab = 'RiskScore',method = 'spearman',
                     marginal.type = 'density')
fig11a
fig11b<-sig_boxplot(dat = tcga.esti.risk[,c("Risk","ImmuneScore")],
                    leg = 'Risk',ylab = 'ImmuneScore',
                    palette = ggsci::pal_simpsons()(9))
fig11b
write.table(tcga.esti.risk,'results/ana/12.module.immnu/tcga.esti.risk.txt',quote = F,row.names = T,sep='\t')
#TIDE
head(tcga.tide)
tcga.tide.esti<-cbind.data.frame(TIDE=tcga.tide[tcga_subtype.cli.risk$Samples,]$TIDE,
                                 RiskScore=tcga_subtype.cli.risk$riskscorez,
                                 Risk=tcga_subtype.cli.risk$Risk)
head(tcga.tide.esti)
rownames(tcga.tide.esti)=tcga_subtype.cli.risk$Samples
fig11c<-sig_boxplot(dat = tcga.tide.esti[,c("Risk","TIDE")],
                    leg = 'Risk',ylab = 'TIDE',
                    palette = ggsci::pal_simpsons()(9))
fig11c
write.table(tcga.tide.esti,'results/ana/12.module.immnu/tcga.tide.esti.txt',quote = F,row.names = T,sep='\t')

#免疫检查点分析
setdiff(c('PDCD1','CD274','CTLA4'),rownames(tcga_tmp_log2_T))

#
icg.genes=c('PDCD1','CD274','CTLA4')
#############"PD-1","PD-L1",'CTLA4'
tcga.icg.genes=t(tcga_tmp_log2_T[icg.genes,])
colnames(tcga.icg.genes)=c("PD-1","PD-L1",'CTLA4')
mg_PlotMutiBoxplot(data=tcga.icg.genes[tcga_subtype.cli.risk$Samples,]
                   ,group = tcga_subtype.cli.risk$Risk
                   ,group_cols = ggsci::pal_simpsons()(9)
                   ,test_method = 'wilcox.test',legend.pos = 't')
colnames(tcga.icg.genes)=c("PD-1**","PD-L1",'CTLA4***')

sam_low=tcga_subtype.cli.risk[which(tcga_subtype.cli.risk$Risk=='Low'),"Samples"]
p1=mg_ridges_plot(data=tcga.icg.genes[sam_low,],xlab = "log2(TPM+1)")+
  theme_bw()+theme(legend.position = 'none')+
  labs(title = "Risk-low")

sam_high=tcga_subtype.cli.risk[which(tcga_subtype.cli.risk$Risk=='High'),"Samples"]

p2=mg_ridges_plot(data=tcga.icg.genes[sam_high,],xlab = "log2(TPM+1)")+
  theme_bw()+theme(legend.position = 'none')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank() )+
  labs(title = "Risk-high")

fig11d=cowplot::plot_grid(p1,p2,align = 'h',
                          nrow = 1,ncol = 2,rel_widths = c(1,1))
fig11d

write.table(tcga.icg.genes,'results/ana/12.module.immnu/tcga.icg.genes.txt',quote = F,row.names = T,sep='\t')

tcga.cytscore=immu_cytScore(tcga_tmp_log2_T)
tcga.cytscore=data.frame(tcga.cytscore,check.names = F)
head(tcga.cytscore)


tcga.cytscore=merge(data.frame(Samples=rownames(tcga.cytscore),
                               score=tcga.cytscore[,1]),
                    tcga_subtype.cli.risk[,c("Samples","Risk","riskscorez")],
                    BY='Samples')
head(tcga.cytscore)
fig11e=sig_boxplot(dat = tcga.cytscore[,c("Risk","score")],
                   leg = 'Risk',
                   palette = ggsci::pal_simpsons()(9),
                   ylab = 'Cytolytic activity')
fig11e

#TRS
TRS.gene=read.delim('origin_datas/TRS.txt',sep='\t',header = F)
TRS.gene=TRS.gene$V1
tcga.trs.ssgsea=ssGSEAScore_by_genes(gene.exp = tcga_tmp_log2_T,genes = as.character(TRS.gene))
tcga.trs.ssgsea=data.frame(t(tcga.trs.ssgsea),check.names = F)
head(tcga.trs.ssgsea)
tcga.trs.ssgsea=merge(data.frame(Samples=rownames(tcga.trs.ssgsea),
                                 score=tcga.trs.ssgsea[,1]),
                      tcga_subtype.cli.risk[,c("Samples","Risk","riskscorez")],
                      BY='Samples')
head(tcga.trs.ssgsea)
write.table(tcga.trs.ssgsea,'results/ana/12.module.immnu/tcga.trs.ssgsea.txt',quote = F,row.names = F,sep='\t')
fig11f=sig_boxplot(dat = tcga.trs.ssgsea[,c("Risk","score")],
                   leg = 'Risk',
                   palette = ggsci::pal_simpsons()(9),
                   ylab = 'TRS score')
fig11f

fig11<-mg_merge_plot(fig11a,fig11b,fig11c,
                     fig11d,fig11e,fig11f,
                     nrow = 2,ncol = 3,labels = LETTERS[1:6])
ggsave('results/ana/12.module.immnu/Fig11.pdf',fig11,height = 9,width = 15)

save.image('project_002.RData')
#补药物
dir.create('results/ana/13.module.drug')
library(pRRophetic)
library(ggplot2)
## Cisplatin,顺铂
set.seed(12345)
predictedPtype_Cisplatin <- pRRopheticPredict(as.matrix(tcga_tmp_log2_T[, rownames(tcga.module.risk$result)])
                                              , "Cisplatin"
                                              , "urogenital_system"
                                              , selection=1
                                              ,dataset = "cgp2016")
predictedPtype_Cisplatin <- data.frame(predictedPtype_Cisplatin)

tcga_durg_ic50_res <- predictedPtype_Cisplatin

drugs <- c("Cisplatin","Erlotinib","Rapamycin","Sunitinib","PHA-665752","MG-132","Paclitaxel","Cyclopamine","AZ628","Sorafenib","VX-680","Imatinib","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine","Z-LLNle-CHO","Dasatinib","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536","BMS-509744","CMK","Pyrimethamine","JW-7-52-1","A-443654","GW843682X","MS-275","Parthenolide","KIN001-135","TGX221","Bortezomib","XMD8-85","Roscovitine","Salubrinal","Lapatinib","Vinorelbine","NSC-87877","QS11","CP466722","Midostaurin","Shikonin","AKT inhibitor VIII","Embelin","Bexarotene","Bleomycin","Phenformin")
dim(tcga_durg_ic50_res)
for (drug in drugs) {
  print(drug)
  set.seed(12345)
  tmpic50 <- pRRopheticPredict(as.matrix(tcga_tmp_log2_T[, rownames(tcga.module.risk$result)])
                               , drug
                               , "urogenital_system"
                               , selection=1
                               , dataset = "cgp2016")
  tmpic50 <- data.frame(tmpic50)
  colnames(tmpic50) <- drug
  tcga_durg_ic50_res <- cbind(tcga_durg_ic50_res, tmpic50)
}

dim(tcga_durg_ic50_res)
colnames(tcga_durg_ic50_res)
colnames(tcga_durg_ic50_res) <- gsub('predictedPtype_', '', colnames(tcga_durg_ic50_res))

dim(tcga_durg_ic50_res)
tcga_durg_ic50_res1 <- cbind.data.frame(tcga_durg_ic50_res[rownames(tcga.module.risk$result),],
                                        riskscore=tcga.module.risk$result$riskscore,
                                        Risk=tcga.module.risk$result$Risk)
head(tcga_durg_ic50_res1)
write.table(tcga_durg_ic50_res1,'results/ana/13.module.drug/tcga_durg_ic50_res.txt',quote = F,sep='\t',row.names = T)

fig15a<-sig_boxplot(dat = tcga_durg_ic50_res1[,c('Risk',drugs[1])],
                    leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[1],')'),
                    palette = ggsci::pal_simpsons()(9)[c(5,2)])
fig15a

fig15b<-sig_boxplot(dat = tcga_durg_ic50_res1[,c('Risk',drugs[2])],
                    leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[2],')'),
                    palette = ggsci::pal_simpsons()(9)[c(5,2)])
fig15b

fig15c<-sig_boxplot(dat = tcga_durg_ic50_res1[,c('Risk',drugs[3])],
                    leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[3],')'),
                    palette = ggsci::pal_simpsons()(9)[c(5,2)])
fig15c

fig15d<-sig_boxplot(dat = tcga_durg_ic50_res1[,c('Risk',drugs[6])],
                    leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[6],')'),
                    palette = ggsci::pal_simpsons()(9)[c(5,2)])
fig15d

fig15e<-sig_boxplot(dat = tcga_durg_ic50_res1[,c('Risk',drugs[9])],
                    leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[9],')'),
                    palette = ggsci::pal_simpsons()(9)[c(5,2)])
fig15e

fig15f<-sig_boxplot(dat = tcga_durg_ic50_res1[,c('Risk',drugs[10])],
                    leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[10],')'),
                    palette = ggsci::pal_simpsons()(9)[c(5,2)])
fig15f

fig15g<-sig_boxplot(dat = tcga_durg_ic50_res1[,c('Risk',drugs[11])],
                    leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[11],')'),
                    palette =ggsci::pal_simpsons()(9)[c(5,2)])
fig15g
fig15h<-sig_boxplot(dat = tcga_durg_ic50_res1[,c('Risk',drugs[15])],
                    leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[15],')'),
                    palette = ggsci::pal_simpsons()(9)[c(5,2)])
fig15h

fig15<-mg_merge_plot(fig15a,fig15b,fig15c,fig15d,
                     fig15e,fig15f,fig15g,fig15h,
                     nrow = 2,ncol = 4)
fig15
ggsave('results/ana/13.module.drug/Fig15.pdf',fig15,height = 7,width = 15)

