library(Seurat)
#load Seurat
setwd("D:/OneDrive - dgist.ac.kr/iNKT/iNKT_2/paper_final")
NCD_seurat <- readRDS("data/NCD_seurat.rds")
NCD_thymus_seurat <- NCD_seurat[,NCD_seurat$anatomy=="Thymus"]
NCD_adipo_seurat <- readRDS("data/NCD_adipo_seurat.rds")
NCD_adipo_seurat <- ScaleData(NCD_adipo_seurat, features = rownames(NCD_adipo_seurat))


library(SingleCellExperiment)
library(KernelKnn)

### Select Highly variable genes (feature selection)
hvg <- VariableFeatures(NCD_adipo_seurat)
hvg <- intersect(hvg, rownames(NCD_seurat))
hvg_mat <- as.matrix(NCD_adipo_seurat@assays$RNA@data)[hvg,]
hvg_test_mat <- as.matrix(NCD_thymus_seurat@assays$RNA@data[hvg,])

indexN <- KernelKnn::knn.index.dist(t(hvg_mat),t(hvg_test_mat),k=5,threads = 4,
                                    method = "pearson_correlation")

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
                   expression= c(rep("NCD_adipo",length(colnames(NCD_adipo_seurat))),as.character(NCD_thymus_seurat$anatomy)))

df_test$expression <- factor(df_test$expression, levels = c("Thymus","NCD_adipo")) 

# Load the "scales" package
require(scales)
# Create vector of default ggplot2 colors
library(ggplot2)
g<-ggplot(df_test,
       aes(x=x, y=y, color = expression)) +
  geom_point(size=0.8) +
  ggtitle("NCD Thymic on NCD adipose umap result")+
  scale_color_manual(values = c("#4472C4","grey79"))+
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

