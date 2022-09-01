library(dplyr)
library(Seurat)
library(patchwork)
# Load the PBMC dataset here I use another one 
ost.data <- Read10X(data.dir ="SRR11955372/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
ost <- CreateSeuratObject(counts =ost.data, project = "ost", min.cells = 3, min.features = 200)
ost


ost[["percent.mt"]] <- PercentageFeatureSet(ost, pattern = "^MT-")
head(ost@meta.data, 5)

VlnPlot(ost, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(ost, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1

plot2 <- FeatureScatter(ost, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2


#计算红细胞比例
HB.genes = c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1", "HBZ")
HB_m <- match(HB.genes,rownames(ost@assays$RNA))
HB.genes = rownames(ost@assays$RNA)[HB_m]
HB.genes = HB.genes[!is.na(HB.genes)]
ost[["percent.HB"]] <- PercentageFeatureSet(ost,features = HB.genes)
#过滤低质量细胞
ost <- subset(ost,subset=nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt<10 & percent.HB<1)



ost <- NormalizeData(ost, normalization.method = "LogNormalize", scale.factor = 10000)


ost <- FindVariableFeatures(ost, selection.method =   
                              "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ost), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ost)
plot2 <- LabelPoints(plot = plot1, points = top10, 
                     repel    = TRUE)
plot1+plot2


all.genes <- rownames(ost)
ost <- ScaleData(ost, features = all.genes)

ost <- RunPCA(ost, features = VariableFeatures(object = ost))

VizDimLoadings(ost, dims = 1:2, reduction = "pca")
DimPlot(ost, reduction = "pca")

DimHeatmap(ost, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(ost, dims = 1:15, cells = 500, balanced = TRUE)

ost <- JackStraw(ost, num.replicate = 100, dim = 50)
ost <- ScoreJackStraw(ost, dims = 1:50)
JackStrawPlot(ost, dims = 1:40)




ElbowPlot(ost)
ElbowPlot(ost, ndims = 50, reduction = "pca")


## 细胞聚类
ost <- FindNeighbors(ost, dims = 1:30)
ost <- FindClusters(ost, resolution = 0.5) 
head(Idents(ost), 5)

# 非线性聚类
ost <- RunUMAP(ost, dims = 1:30)
head(ost@reductions$umap@cell.embeddings)  #提取UMAP坐标值
DimPlot(ost, reduction = "umap")  #降维方法首选umap，再tsne，最后再pca

cluster2.markers <- FindMarkers(ost, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

cluster5.markers <- FindMarkers(ost, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

ost.markers <- FindAllMarkers(ost, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ost.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(ost, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(ost, features = c("COL2A1", "S100A1"))
VlnPlot(ost, features = c("CXCL14", "SFRP2"), slot = "counts", log = TRUE)

FeaturePlot(ost, features = c("PTN", "COL2A1", "SFRP2", "HLA-DRA", "ASPM", "CAV1", "SPARCL1","TAGLN"))

library(celldex)
library(SingleR)
hpca.ref <- HumanPrimaryCellAtlasData()
pred <- SingleR(test = ost@assays$RNA@data, ref = hpca.ref, labels = hpca.ref$label.main,  clusters = ost@active.ident)

table(pred$labels)
pdf("cell_type.pdf")
plotScoreHeatmap(pred)

plotDeltaDistribution(pred, ncol = 3)
new.cluster.ids <- pred$labels
names(new.cluster.ids) <- levels(ost)

ost <- RenameIdents(ost, new.cluster.ids)

DimPlot(ost, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

new.cluster.ids
