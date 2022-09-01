library(dplyr)
library(Seurat)
library(patchwork)

dir0 <- list.files('./',pattern='SRR')
scRNA.list <- list()
dir <-list()
for(i in 1:length(dir0)){
  dir[[i]]<-list.dirs(paste('./',dir0[i],'/', sep = ""))
  counts=Read10X(data.dir = dir[[i]][-1])
  scRNA.list[[i]]<- CreateSeuratObject(counts,project=dir0[i],min.cells = 3,min.features = 200)
}



for(i in 1:length(dir)){
  scRNA.list[[i]]<- RenameCells(scRNA.list[[i]],add.cell.id = dir0[i])
  scRNA.list[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNA.list[[i]], pattern = "^MT-")
  scRNA.list[[i]] <- subset(scRNA.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
  scRNA.list[[i]]<- NormalizeData(scRNA.list[[i]])
  scRNA.list[[i]]<- FindVariableFeatures(scRNA.list[[i]],selection.method='vst')
  
}

#锚点法
#找到不同细胞共有的高变基因
features <- SelectIntegrationFeatures(object.list=scRNA.list)
#找到锚定点
scRNA.anchors <- FindIntegrationAnchors(object.list=scRNA.list,anchor.features=features)
#进行数据整合
scRNA1 <- IntegrateData(anchorset = scRNA.anchors)
#接下来进行ScaleData、PCA、TSNE、UMAP等处理（略）

all.genes <- rownames(scRNA1)
scRNA1 <- ScaleData(scRNA1, features = all.genes)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method =   
                                     "vst", nfeatures = 2000)
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(object = scRNA1))

## 细胞聚类
scRNA1 <- FindNeighbors(scRNA1, dims = 1:30)
scRNA1 <- FindClusters(scRNA1, resolution = 0.5) 
head(Idents(scRNA1), 5)


scRNA1 <- RunUMAP(scRNA1, dims = 1:30)
head(scRNA1@reductions$umap@cell.embeddings)  #提取UMAP坐标值
DimPlot(scRNA1, reduction = "umap")  #降维方法首选umap，再tsne，最后再pca

DimPlot(scRNA1, reduction = "umap",group.by = 'orig.ident')
DimPlot(scRNA1, reduction = "umap",split.by = 'orig.ident')

#harmony
scRNA_harmony <- merge(scRNA.list[[1]], y = scRNA.list[2:length(scRNA.list)])
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
#group.by.vars指定分组依据
scRNA_harmony <- RunHarmony(scRNA_harmony,group.by.vars = "orig.ident")

scRNA_harmony <- RunUMAP(scRNA_harmony,reduction = "harmony",dims=1:40)
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction = "harmony",dims=1:40) %>% FindClusters(resolution = 0.5)

plot1=DimPlot(scRNA_harmony,reduction = "umap",label=T)
plot1



Merge_data <- FindNeighbors(Merge_data, dims = 1:50)
Merge_data <- FindClusters(Merge_data, resolution = 0.5) 
head(Idents(Merge_data), 5)

Merge_data <- RunUMAP(Merge_data, dims = 1:50)
head(Merge_data@reductions$umap@cell.embeddings)  #提取UMAP坐标值
DimPlot(Merge_data, reduction = "umap")  #降维方法首选umap，再tsne，最后再pca



#此处省略过滤线粒体与红细胞的质控过程
#第三步 进行归一化处理并寻找高变基因

# Load the ost dataset here I use another one 
ost.data <- Read10X(data.dir ="SRR11955372/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
ost <- CreateSeuratObject(counts =ost.data, project = "ost", min.cells = 3, min.features = 200)
ost


ost[["percent.mt"]] <- PercentageFeatureSet(ost, pattern = "^MT-")

VlnPlot(ost, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(ost, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1

plot2 <- FeatureScatter(ost, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

ost <- NormalizeData(ost, normalization.method = "LogNormalize", scale.factor = 10000)


ost <- FindVariableFeatures(ost, selection.method =   
                              "vst", nfeatures = 2000)

library(devtools)
install_github("immunogenomics/harmony")
library(harmony)
