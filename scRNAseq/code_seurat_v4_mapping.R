.libPaths()

library(Seurat)
seurat <- readRDS("seurat.rds")
ncd_adipose_seurat <- readRDS("NCD_adipo_seurat.rds")
NCD_seurat <- readRDS("NCD_seurat.rds")
seurat_V4 <- seurat_subset_recluster(seurat,cells = colnames(seurat),
                                     PCA = 20,n=1000)

tissue.list <- SplitObject(seurat_V4, split.by = "tissue")
for (i in 1:length(tissue.list)) {
  tissue.list[[i]] <- NormalizeData(tissue.list[[i]], verbose = FALSE)
  tissue.list[[i]] <- FindVariableFeatures(tissue.list[[i]], selection.method = "vst", nfeatures = 1000,
                                           verbose = FALSE)
}
reference.list <- tissue.list[c("Adipose","Thymus")]
tissue.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)
tissue.integrated <- IntegrateData(anchorset = tissue.anchors, dims = 1:20)
saveRDS(tissue.integrated,"tissue.integrated_V4.rds")

reference.list <- tissue.list[c("Adipose","Thymus","Spleen","Liver","LymphNode")]
tissue.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)
tissue.integrated <- IntegrateData(anchorset = tissue.anchors, dims = 1:20)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(tissue.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
tissue.integrated <- ScaleData(tissue.integrated, verbose = FALSE)
tissue.integrated <- RunPCA(tissue.integrated, npcs = 20, verbose = FALSE)
tissue.integrated <- RunUMAP(tissue.integrated, reduction = "pca", dims = 1:20)
DimPlot_eunseo(tissue.integrated,reduction = "umap",
               split.by = "tissue")
DimPlot_eunseo(tissue.integrated,reduction = "umap",
               cells.highlight = colnames(ncd_adipose_seurat)[ncd_adipose_seurat$final_annot=="As_NKT1"])
saveRDS(tissue.integrated,"data/alltissue.integrated_V4.rds")
library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(tissue.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
tissue.integrated <- ScaleData(tissue.integrated, verbose = FALSE)
tissue.integrated <- RunPCA(tissue.integrated, npcs = 20, verbose = FALSE)
tissue.integrated <- RunUMAP(tissue.integrated, reduction = "pca", dims = 1:20, verbose = FALSE)
p1 <- DimPlot(tissue.integrated, reduction = "umap", group.by = "tissue")
p2 <- DimPlot(tissue.integrated, reduction = "umap", group.by = "annot", label = TRUE, repel = TRUE) +
  NoLegend()
p1 + p2

all_df <- data.frame()

for(i in names(tissue.list)[1:4]){
  tissue.query <- tissue.list[[i]]
  tissue.anchors <- FindTransferAnchors(reference = ncd_adipose_seurat, query = tissue.query,
                                        dims = 1:30, reference.reduction = "pca")
  predictions <- TransferData(anchorset = tissue.anchors, refdata = ncd_adipose_seurat$final_annot,
                              dims = 1:30)
  tissue.query <- AddMetaData(tissue.query, metadata = predictions)
  
  tmp_df <- as.data.frame(tissue.query$predicted.id)
  tmp_df$tissue <- i
  all_df <- rbind(tmp_df,all_df)
}

library(dplyr)
colnames(all_df)[1] <- "predicted_celltype"
all_df$predicted_celltype <- as.factor(all_df$predicted_celltype)
percent_df <- all_df %>% group_by(tissue,predicted_celltype,.drop = F) %>% summarise(count=n())  %>% mutate(perc=count/sum(count))
write.csv(as.data.frame(percent_df),"data/all_tissue_mapping_result.csv")
library(ggplot2)
g <- ggplot(percent_df[,], aes(x=tissue, y=perc*100, fill=predicted_celltype)) + 
  geom_bar(stat = "identity",position = "dodge") +
  scale_fill_manual(values = celltype_colors)+
  theme(  # Remove panel border 
    panel.border = element_blank(),  
    # Remove panel grid linestkfdl
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=10),
    axis.text.x = element_text(size = 10,angle = 90,colour = "black"))
g


tissue.query <- tissue.list[["Thymus"]]
tissue.anchors <- FindTransferAnchors(reference = ncd_adipose_seurat, query = tissue.query,
                                      dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = tissue.anchors, refdata = ncd_adipose_seurat$final_annot,
                            dims = 1:30)
tissue.query <- AddMetaData(tissue.query, metadata = predictions)
thymus_query <- tissue.query

ncd_adipose_seurat_V4 <- seurat_subset_recluster(ncd_adipose_seurat,cells = colnames(ncd_adipose_seurat),
                                     PCA = 30,n=1000)
DimPlot_eunseo(ncd_adipose_seurat_V4,
               group.by = "final_annot")
ncd_adipose_seurat_V4 <- RunUMAP(ncd_adipose_seurat_V4,return.model = T,
                                 dims = 1:30)
check.query <- MapQuery(anchorset = tissue.anchors, reference = ncd_adipose_seurat_V4, query = tissue.query,
                           refdata = list(celltype = "annot"), reference.reduction = "pca", reduction.model = "umap")
DimPlot(check.query, reduction = "ref.umap",
        group.by = "predicted.id")
DimPlot_eunseo(ncd_adipose_seurat_V4,group.by = "final_annot")
DimPlot_eunseo(check.query)
thymus_query$annot <- as.factor(thymus_query$annot)
thymus_query$predicted.id <- as.factor(thymus_query$predicted.id)
thymus_df <- thymus_query@meta.data %>% dplyr::group_by(annot,predicted.id,.drop = FALSE) %>% 
  summarise(count=n()) %>% mutate(perc=count/sum(count))

library(ggplot2)
g <-ggplot(thymus_df[,], aes(x=annot, y=perc*100, fill=predicted.id)) +
  geom_bar(stat = "identity",position = "dodge") +
  #scale_fill_manual(values = cols_df$colors)+
  scale_fill_manual(values = celltype_colors)+
  theme(  # Remove panel border 
    panel.border = element_blank(),  
    # Remove panel grid linestkfdl
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=10),
    axis.text.x = element_text(size = 10,angle = 90,colour = "black"))
g


