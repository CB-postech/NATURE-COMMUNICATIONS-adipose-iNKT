.libPaths(Sys.getenv("iNKT_R3"))
.libPaths()
library(SingleCellExperiment)
library(DropletUtils)
library(scater)
iNKT_R3
setwd("D:/OneDrive - dgist.ac.kr/iNKT/iNKT_2/raw_data")
sample_info_list <- list.files("D:/OneDrive - dgist.ac.kr/iNKT/iNKT_2/raw_data/")
plotdir <-"D:/OneDrive - dgist.ac.kr/iNKT/iNKT_2/raw_data/"
sample_info_list <-sample_info_list[1:5]

#DropletUtils
for( i in sample_info_list){
  dir <- paste0("D:/OneDrive - dgist.ac.kr/iNKT/iNKT_2/raw_data/",i,"/raw_feature_bc_matrix")
  sce <- read10xCounts(dir)
  rownames(sce) = uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
  colnames(sce) = sce$Barcode
  
  my.counts=counts(sce)
  br.out <- barcodeRanks(my.counts)
  png(paste0(i,".png"))
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  dev.off()
  set.seed(100)
  e.out <- emptyDrops(my.counts)
  is.cell <- e.out$FDR <= 0.05
  print(paste0("The number of QC positive cells of ",i," is ",sum(is.cell, na.rm=TRUE)))
  is.cell[is.na(is.cell)] <- FALSE
  
  names(is.cell) <- colnames(sce)
  sce$cells_kept <- is.cell
  
  sce <- sce[, sce$cells_kept == T]
  sce
  assign(paste0(i,"_sce"),sce)
}
library(RColorBrewer)
####Set QC filter : log(counts)>2.5 mito <10
for(i in sample_info_list){
  
  objectname = paste0(i,"_sce")
  sce = eval(parse(text=objectname))
  sce1 <- qc_function(sce, objectname = objectname,
                           total_counts_cutoff = 2.5, mt_percent_cutoff = 10)
  sce.qc <- qc_function(sce, objectname = objectname,
                             total_counts_cutoff = 2.5, mt_percent_cutoff = 10, save = TRUE)
  assign(paste0(i,"_sce"), sce1)
  assign(paste0(i,"_sce", ".qc"), sce.qc)
}
for(i in sample_info_list){
  objectname = paste0(i,"_sce.qc")
  sce = eval(parse(text=objectname))
  saveRDS(sce,paste0(i,"sce_qc_count2_5_mito10.rds"))
}

####Set QC filter : log(counts)>3 mito <10
for(i in sample_info_list){
  
  objectname = paste0(i,"_sce")
  sce = eval(parse(text=objectname))
  sce1 <- qc_function(sce, objectname = objectname,
                           total_counts_cutoff = 3, mt_percent_cutoff = 10)
  sce.qc <- qc_function(sce, objectname = objectname,
                             total_counts_cutoff = 3, mt_percent_cutoff = 10, save = TRUE)
  assign(paste0(i,"_sce"), sce1)
  assign(paste0(i,"_sce", ".qc"), sce.qc)
}

for(i in sample_info_list){
  objectname = paste0(i,"_sce.qc")
  sce = eval(parse(text=objectname))
  saveRDS(sce,paste0(i,"sce_qc_count3_mito10.rds"))
}

###Save Droplet QC SCE
for(i in sample_info_list){
  objectname = paste0(i,"_sce")
  sce = eval(parse(text=objectname))
  saveRDS(sce,paste0(i,"sce_droplet_QC.rds"))
}
