# Metabolic pathway scoring for Seurat objects using scMetabolism.
#
# This script computes metabolism scores (KEGG / REACTOME) for a selected subset of
# T cell clusters from an integrated Seurat object.
#
# Inputs
# - A saved Seurat object (`.Rdata`) with a `celltype` column in metadata.
#
# Outputs
# - Tab-delimited score matrices written to the specified results folders.
#
# Notes
# - The analysis uses `sc.metabolism.Seurat(..., method = "VISION")` which may require
#   substantial CPU/memory; adjust `ncores` appropriately.

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
#load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/seurat.integrated.T_refined_cluster.Rdata")
#seurat.integrated.T_cell_sub=seurat.integrated.T_cell_sub[,which(seurat.integrated.T_cell_sub@meta.data$celltype %in% sort(unique(seurat.integrated.T_cell_sub@meta.data$celltype))[12:19])]

#seurat.integrated.T_cell_sub<-sc.metabolism.Seurat(obj = seurat.integrated.T_cell_sub, method = "VISION", imputation = F, ncores = 20, metabolism.type = "KEGG")
#seurat.integrated.T_cell_sub@assays$METABOLISM$score

#makedirs("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/METABOLISM_KEGG_CD8T/")
#write.table(seurat.integrated.T_cell_sub@assays$METABOLISM$score,file="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/METABOLISM_KEGG_CD8T/METABOLISM_score_KEGG.xls",row.names = T,col.names = T,sep="\t")


#load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/seurat.integrated.T_refined_cluster.Rdata")
#seurat.integrated.T_cell_sub=seurat.integrated.T_cell_sub[,which(seurat.integrated.T_cell_sub@meta.data$celltype %in% sort(unique(seurat.integrated.T_cell_sub@meta.data$celltype))[12:19])]
#seurat.integrated.T_cell_sub<-sc.metabolism.Seurat(obj = seurat.integrated.T_cell_sub, method = "VISION", imputation = F, ncores = 20, metabolism.type = "REACTOME")
#seurat.integrated.T_cell_sub@assays$METABOLISM$score

#makedirs("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/METABOLISM_REACTOME_CD8T/")
#write.table(seurat.integrated.T_cell_sub@assays$METABOLISM$score,file="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/METABOLISM_REACTOME_CD8T/METABOLISM_score_REACTOME.xls",row.names = T,col.names = T,sep="\t")





# -----------------------------
# KEGG metabolism scoring
# -----------------------------
load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/seurat.integrated.T_refined_cluster.Rdata")
# Example subset: restrict to a set of T cell subtypes and a TLS stratum.
seurat.integrated.T_cell_sub=seurat.integrated.T_cell_sub[,which(seurat.integrated.T_cell_sub@meta.data$celltype %in% sort(unique(seurat.integrated.T_cell_sub@meta.data$celltype))[9:16] & seurat.integrated.T_cell_sub@meta.data$TLS=="With_TLS")]
seurat.integrated.T_cell_sub<-sc.metabolism.Seurat(obj = seurat.integrated.T_cell_sub, method = "VISION", imputation = F, ncores = 20, metabolism.type = "KEGG")
#seurat.integrated.T_cell_sub@assays$METABOLISM$score

# Save KEGG score matrix (cells x pathways) for downstream plotting/statistics.
makedirs("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/METABOLISM_KEGG_CD8T/")
write.table(seurat.integrated.T_cell_sub@assays$METABOLISM$score,file="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/METABOLISM_KEGG_CD8T/METABOLISM_score_KEGG_WithTLS.xls",row.names = T,col.names = T,sep="\t")


# -----------------------------
# REACTOME metabolism scoring
# -----------------------------
load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/seurat.integrated.T_refined_cluster.Rdata")
seurat.integrated.T_cell_sub=seurat.integrated.T_cell_sub[,which(seurat.integrated.T_cell_sub@meta.data$celltype %in% sort(unique(seurat.integrated.T_cell_sub@meta.data$celltype))[9:16] & seurat.integrated.T_cell_sub@meta.data$TLS=="With_TLS")]
seurat.integrated.T_cell_sub<-sc.metabolism.Seurat(obj = seurat.integrated.T_cell_sub, method = "VISION", imputation = F, ncores = 20, metabolism.type = "REACTOME")
#seurat.integrated.T_cell_sub@assays$METABOLISM$score

# Save REACTOME score matrix.
makedirs("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/METABOLISM_REACTOME_CD8T/")
write.table(seurat.integrated.T_cell_sub@assays$METABOLISM$score,file="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/METABOLISM_REACTOME_CD8T/METABOLISM_score_REACTOME_WithTLS.xls",row.names = T,col.names = T,sep="\t")
