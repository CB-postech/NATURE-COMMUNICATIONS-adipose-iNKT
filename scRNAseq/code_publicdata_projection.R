library(Seurat)
library(SingleCellExperiment)
library(scran)
library(scater)
library(DropletUtils)

setwd("public_data/")
#load Seurat
sce_merge <- cbind(GSM4908722_NKT_mm_Liver_sce,GSM4908723_NKT_mm_Liver_sce,
                   GSM4908724_NKT_mm_LymphNode_sce,GSM4908725_NKT_mm_LymphNode_sce,
                   GSM4908726_NKT_mm_Spleen_sce,GSM4908727_NKT_mm_Spleen_sce)
keep_feature <- rowSums(counts(sce_merge) != 0) > 0
sce_merge <- sce_merge[keep_feature, ]
### filtering genes
library(scran)
set.seed(123)
clusters <- quickCluster(sce_merge, method="igraph")
sce_merge <- computeSumFactors(sce_merge, cluster=clusters)
sce_merge <- logNormCounts(sce_merge)
dec <- modelGeneVar(sce_merge)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

hvg <- getTopHVGs(dec,n=1000)
#hvg <- getTopHVGs(dec,fdr.threshold = 0.05)
length(hvg)

library(Seurat)
seurat <- as.Seurat(sce_merge)
VariableFeatures(seurat) <- hvg
seurat <- ScaleData(seurat)
seurat@reductions$PCA_coldata <- NULL
plot(seurat@reductions$pca@stdev)

PCA = 20
seurat<- RunPCA(seurat, npcs = PCA)
seurat<- FindNeighbors(seurat, reduction="pca", dims= 1:PCA)
seurat<- FindClusters(seurat)
seurat<- RunUMAP(seurat, dims = 1:PCA, seed.use = 42)
seurat$sample <- substring(colnames(seurat),18)
DimPlot_eunseo(seurat,group.by = "sample")
seurat$tissue <- "Lymphnode"
seurat$tissue[grep("Liver",seurat$sample)]<- "Liver"
seurat$tissue[grep("Spleen",seurat$sample)]<- "Spleen"
DimPlot_eunseo(seurat,group.by = "tissue")
spleen_seurat <- seurat_subset_recluster(seurat,cells = colnames(seurat)[seurat$tissue=="Spleen"],n=1000)
liver_seurat <- seurat_subset_recluster(seurat,cells = colnames(seurat)[seurat$tissue=="Liver"],n=1000)
lymphnode_seurat <- seurat_subset_recluster(seurat,cells = colnames(seurat)[seurat$tissue=="Lymphnode"],n=1000)

library(SingleCellExperiment)
library(KernelKnn)
seurat <- spleen_seurat
### Select Highly variable genes (feature selection)
NCD_adipo_seurat <- ScaleData(NCD_adipo_seurat, features = rownames(NCD_adipo_seurat))
hvg <- VariableFeatures(NCD_adipo_seurat)
hvg <- intersect(hvg, rownames(seurat))
hvg_mat <- as.matrix(NCD_adipo_seurat@assays$RNA@data)[hvg,]
hvg_test_mat <- as.matrix(seurat@assays$RNA@data[hvg,])
# 
# hvg_mat <- as.matrix(NCD_adipo_seurat@assays$RNA@scale.data)[hvg,]
# hvg_test_mat <- as.matrix(public_all_seurat@assays$RNA@scale.data[hvg,])

indexN <- KernelKnn::knn.index.dist(t(hvg_mat),t(hvg_test_mat),k=5,threads = 4,
                                    method = "pearson_correlation")

# saveRDS(indexN,"NCD_all_on_NCD_adipo_default_parameter.rds")
# saveRDS(indexN,"NCD_all_on_NCD_adipo_pearson_parameter.rds")

# indexN <- readRDS("I:/SNU_adipo_infla/NCD_all_on_NCD_adipo_pearson_parameter.rds")
iN2 <- indexN$test_knn_idx
rownames(iN2) <- colnames(hvg_test_mat)

NCD_adipo_umap <- as.matrix(NCD_adipo_seurat@reductions$umap@cell.embeddings)

iN3 <- apply(iN2, 2, function(x) colnames(hvg_mat)[x])
idumapx <- apply(iN3, 2, function(x) NCD_adipo_umap[x,1])
idumapy <- apply(iN3, 2, function(x) NCD_adipo_umap[x,2])

prjumapx <- rowMeans(idumapx)
prjumapy <- rowMeans(idumapy)
names(prjumapx) <- rownames(iN2)
names(prjumapy) <- rownames(iN2)

df_test=data.frame(x=rbind(as.matrix(NCD_adipo_umap[,1]),as.matrix(prjumapx)), 
                   y=rbind(as.matrix(NCD_adipo_umap[,2]),as.matrix(prjumapy)), 
                   expression= c(rep("NCD_adipo",length(colnames(NCD_adipo_seurat))),rep("public",ncol(seurat))))
df_test$expression <- factor(df_test$expression, levels=c("NCD_adipo","public")) 

# Load the "scales" package
require(scales)
# Create vector of default ggplot2 colors
library(ggplot2)
df_test
dev.off()
g<-ggplot(df_test,
       aes(x=x, y=y, color = expression)) +
  geom_point(size=0.8) +
  ggtitle("Public data on NCD adipose umap result")+
  scale_color_manual(values = c("grey79","#4DB34F"))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  #scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  #scale_colour_gradient(low = "grey", high = "red") +
  #scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
  #breaks=c(min(celltype_data),max(celltype_data)),
  #labels=c(round(as.numeric(min(celltype_data)),3),round(as.numeric(max(celltype_data)), digits = 2))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())
g
spleen_color <- "#4DB34F"
lymph_node_color <- "#9B7CC9"
liver_color <- "#AF784F"
