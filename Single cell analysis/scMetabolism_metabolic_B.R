# Metabolic pathway scoring for B cell Seurat objects using scMetabolism.
#
# This script loads an integrated B cell Seurat object, computes KEGG and REACTOME
# metabolism scores using the VISION method, and writes the score matrices to disk.
#
# Inputs
# - `/NFS_home/.../seurat.integrated.B_refined_cluster.Rdata`
#
# Outputs
# - `METABOLISM_score_KEGG.xls`
# - `METABOLISM_score_REACTOME.xls`
#
# Notes
# - `sc.metabolism.Seurat` adds a `METABOLISM` assay containing the score matrix.
# - Adjust `ncores` if running on a smaller machine.

library(ggplot2)
library(Seurat)
library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(6,"Dark2")
colorl <- rep(color,each=10)

# Utility: create directories if they do not exist.
makedirs <- function(dir_list){ 
  for (i_dir in dir_list){ 
    if (!dir.exists(i_dir) ){ 
      dir.create(i_dir, recursive = TRUE) 
    } 
  } 
}
library(scMetabolism)
library(ggplot2)
library(rsvd)
# -----------------------------
# KEGG metabolism scoring
# -----------------------------
load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/seurat.integrated.B_refined_cluster.Rdata")
#seurat.integrated.T_cell_sub=seurat.integrated.T_cell_sub[,rownames(seurat.integrated.T_cell_sub@meta.data[which(! seurat.integrated.T_cell_sub@meta.data$celltype %in% c("B06_B-APOC1","B08_B-CXCL8")),])]
seurat.integrated.T_cell_sub<-sc.metabolism.Seurat(obj = seurat.integrated.T_cell_sub, method = "VISION", imputation = F, ncores = 20, metabolism.type = "KEGG")
#seurat.integrated.T_cell_sub@assays$METABOLISM$score

# Save KEGG score matrix.
makedirs("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/METABOLISM_KEGG/")
write.table(seurat.integrated.T_cell_sub@assays$METABOLISM$score,file="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/METABOLISM_KEGG/METABOLISM_score_KEGG.xls",row.names = T,col.names = T,sep="\t")


# -----------------------------
# REACTOME metabolism scoring
# -----------------------------
load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/seurat.integrated.B_refined_cluster.Rdata")
seurat.integrated.T_cell_sub<-sc.metabolism.Seurat(obj = seurat.integrated.T_cell_sub, method = "VISION", imputation = F, ncores = 20, metabolism.type = "REACTOME")
#seurat.integrated.T_cell_sub@assays$METABOLISM$score

# Save REACTOME score matrix.
makedirs("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/METABOLISM_REACTOME/")
write.table(seurat.integrated.T_cell_sub@assays$METABOLISM$score,file="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/METABOLISM_REACTOME/METABOLISM_score_REACTOME.xls",row.names = T,col.names = T,sep="\t")









