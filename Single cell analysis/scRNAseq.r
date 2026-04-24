# Single-cell RNA/BCR/TCR analysis workflow (Seurat + DoubletFinder + Monocle).
#
# This script is written as a sequential notebook-style pipeline (not as functions).
# The intent is to be run top-to-bottom in an R session on an HPC filesystem.
#
# High-level stages
# 1) Load CellRanger count matrices and create one Seurat object per sample
# 2) Merge samples into a single Seurat object and add metadata (sample/barcode)
# 3) Doublet removal with DoubletFinder (per sample)
# 4) Normalization, integration, scaling, dimensionality reduction, clustering
# 5) Cell-type annotation and visualization
# 6) Subsetting (e.g., B cell re-clustering) and marker analysis/heatmaps
# 7) Trajectory inference with Monocle (DDRTree + ordering)
#
# Inputs
# - CellRanger output folders under:
#   `/NFS_home/NFS_home_6/zhouchi/scRNA/cellranger_out/count/<Txxxxxx>/outs/filtered_feature_bc_matrix/`
#
# Outputs (examples; many paths are hard-coded below)
# - Integrated Seurat objects saved via `saveRDS` / `save`
# - UMAP/TSNE plots (`ggsave`)
# - Heatmaps and Monocle objects

suppressMessages(library(Seurat))
suppressMessages(library(DoubletFinder))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
#suppressMessages(library(clusterProfiler))
library(msigdbr)
#library(velocyto.R)
library(DoubletFinder)

# -------------------------------
# 1) Load single-cell matrices
# -------------------------------
setwd('/NFS_home/NFS_home_6/zhouchi/scRNA/cellranger_out/count')
dirs <- list.files(path = '.', pattern = '^T[0-9-]+$')
#dirs <- c('T220799','T220826','T220831','T221053')
for(x in dirs){
    name <- x
    pathway=paste0(x,'/outs/filtered_feature_bc_matrix/')
    # Read sparse matrix (genes x cells) from CellRanger output.
    cts <- ReadMtx(mtx = paste0(pathway,'matrix.mtx.gz'),
                   features = paste0(pathway,'features.tsv.gz'),
                   cells = paste0(pathway,'barcodes.tsv.gz'))
    # Create a Seurat object and assign it to a variable with the same name as the sample ID.
    assign(name, CreateSeuratObject(counts=cts))
}
merged.Seurat <- merge(T213082, y = c(T213087,T220002-3,T220008,T220010,T220011,T220012-4,T220421,T220442,T220452,T220697,T220719,T220721,T220723,T220799,T220826,T220831,T221053),
      add.cell.ids = dirs,
      project='ICC')

merged.Seurat$cell <- row.names(merged.Seurat@meta.data)
merged.Seurat@meta.data$sample=substr(merged.Seurat@meta.data$cell,start = 1,stop = (str_length(merged.Seurat@meta.data$cell)-19))
merged.Seurat@meta.data$barcode=substr(merged.Seurat@meta.data$cell,start = (str_length(merged.Seurat@meta.data$cell)-17),stop = str_length(merged.Seurat@meta.data$cell))


# --------------------------------------------
# 2) Doublet detection and removal (per sample)
# --------------------------------------------
# This loop performs DoubletFinder on each sample subset, then merges results back into
# the combined object. Parameters (pK, expected doublet rate, etc.) are set below in the
# downstream code blocks.
for (s in dirs){
    seurat_sample=merged.Seurat[,which(merged.Seurat@meta.data$sample %in% s)]
    x=seurat_sample
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x <- ScaleData(x, features = rownames(x))
    x <- RunPCA(x, features = VariableFeatures(x),npcs = 50)
    x <- FindNeighbors(x, dims = 1:30)
    x <- FindClusters(x, resolution = 0.5)
    x <- RunUMAP(x, dims = 1:30)
    #x <- RunTSNE(x, dims = 1:20)
    sweep.res.list <- paramSweep(x, PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    DoubletRate= ncol(x)*16*1e-6 
    homotypic.prop <- modelHomotypic(x$seurat_clusters)
    nExp_poi <- round(DoubletRate*ncol(x))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    x <- doubletFinder(x, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
    colnames(x@meta.data)[ncol(x@meta.data)]="DoubletFinder"
    DimPlot(x,reduction = "umap",pt.size = 2,group.by = "DoubletFinder")
    seurat_sample=seurat_sample[,which(x@meta.data$DoubletFinder=='Singlet')]
    save(seurat_sample, file = paste0("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/Singlet/",s,"_remove_doublet.RData"))
}


# Merge and quality control

sample_list=dirs
load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/Singlet/T213082_remove_doublet.RData")
merged.Seurat=seurat_sample

for (s in sample_list[2:17]){
    load(paste0("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/Singlet/",s,"_remove_doublet.RData"))
    Seurat=seurat_sample
    merged.Seurat=merge(merged.Seurat,Seurat)
}

seurat_mtx <- merged.Seurat
seurat_mtx$log10GenesPerUMI <-log10(seurat_mtx$nFeature_RNA) / log10(seurat_mtx$nCount_RNA)
seurat_mtx[['percent.mt']] <- PercentageFeatureSet(seurat_mtx, pattern="^MT-")
seurat_mtx[['percent.ribo']] <- PercentageFeatureSet(seurat_mtx, pattern="^RP[SL]")
rb.genes <- rownames(seurat_mtx)[grep("^RP[SL]",rownames(seurat_mtx))]
#MT.genes <- rownames(seurat_mtx)[grep("^MT",rownames(seurat_mtx))]
#HSP.genes <- c(rownames(seurat_mtx)[grep("^HSP",rownames(seurat_mtx))],rownames(seurat_mtx)[grep("^DNAJ",rownames(seurat_mtx))])
remove_gene=c(rb.genes)
gene_unexclude=setdiff(rownames(seurat_mtx),remove_gene)
seurat_mtx=seurat_mtx[gene_unexclude,]
seurat_mtx_filtered <- subset(seurat_mtx, subset = nFeature_RNA > 100 & nCount_RNA > 500  & percent.mt < 15)
seurat_mtx_filtered
save(seurat_mtx_filtered, file = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_bigcluster/seurat_mtx_filtered_before_batch.RData")


# remove batch effect and clustering
load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_bigcluster/seurat_mtx_filtered_before_batch.RData")
obj.list <- SplitObject(seurat_mtx_filtered, split.by = 'sample')
for(i in 1:length(obj.list)){
    obj.list[[i]] <- NormalizeData(object=obj.list[[i]])
    obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}
features <- SelectIntegrationFeatures(object.list=obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                 anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors)
seurat.integrated <- ScaleData(seurat.integrated)
seurat.integrated <- RunPCA(seurat.integrated)
seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30)
seurat.integrated <- FindClusters(seurat.integrated, resolution = c(0.6,0.8,1.0))
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:30)
seurat.integrated <- RunTSNE(seurat.integrated, dims = 1:30)
# save results
save(seurat.integrated, file = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_bigcluster/seurat.integrated.RData")


# Major celltype definition
B_cell=c(12,13,24,33)
T_cell=c(0,1,3,4,5,9,11,21,22,23,26,29,30,32)
NK_cell=c(15,16,18,25)
Myeloid=c(2,6,7,8,10,14,17,20,27,28)
Epithelial=c(19,31)

current.cluster.ids <- c(B_cell,
                         T_cell,
                         NK_cell,
                         Myeloid,
                         Epithelial)

new.cluster.ids <- c(rep("B_cell",length(B_cell)),
                     rep("T_cell",length(T_cell)),
                     rep("NK_cell",length(NK_cell)),
                     rep("Myeloid",length(Myeloid)),
                     rep("Epithelial",length(Epithelial)))


seurat.integrated@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(seurat.integrated@meta.data$`integrated_snn_res.1`)), from = current.cluster.ids, to = new.cluster.ids)
plotCB=as.data.frame(seurat.integrated@meta.data)[,"cell"]
DimPlot(seurat.integrated, reduction = "umap", group.by = "celltype", pt.size=0.5,cells = plotCB,label=TRUE,raster=FALSE)
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_bigcluster/umap_big_cluster.pdf", width = 15, height = 12, units = 'cm')
saveRDS(seurat.integrated,file = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_bigcluster/seurat.integrated.bigCluster.rds") #Save the Seurat object so it can be loaded directly next time.


# B cell recluster and celltype annotation
#load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/seurat.integrated_B_filter.RData")

seurat.integrated=seurat.integrated[,which(seurat.integrated@meta.data$`integrated_snn_res.1`=="12" | seurat.integrated@meta.data$`integrated_snn_res.1`=="13" | seurat.integrated@meta.data$`integrated_snn_res.1`=="24" | seurat.integrated@meta.data$`integrated_snn_res.1`=="33")]

scRNA_T_cell <- seurat.integrated
scRNA_T_cell$sample=gsub("T220002-3","T213082",scRNA_T_cell$sample)
scRNA_T_cell$sample=gsub("T220008","T213082",scRNA_T_cell$sample)
scRNA_T_cell$sample=gsub("T220010","T213082",scRNA_T_cell$sample)
scRNA_T_cell$sample=gsub("T220421","T213082",scRNA_T_cell$sample)
scRNA_T_cell$sample=gsub("T213087","T213082",scRNA_T_cell$sample)
scRNA_T_cell$sample=gsub("T220719","T213082",scRNA_T_cell$sample)
scRNA_T_cell$sample=gsub("T220826","T213082",scRNA_T_cell$sample)
scRNA_T_cell_sub=CreateSeuratObject(counts = scRNA_T_cell@assays$RNA@counts,
                       meta.data = scRNA_T_cell@meta.data)

obj.list <- SplitObject(scRNA_T_cell_sub, split.by = 'sample')
for(i in 1:length(obj.list)){
    obj.list[[i]] <- NormalizeData(object=obj.list[[i]])
    obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}
features <- SelectIntegrationFeatures(object.list=obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                 anchor.features = features)
seurat.integrated.T_cell_sub <- IntegrateData(anchorset = anchors)

seurat.integrated.T_cell_sub <- ScaleData(seurat.integrated.T_cell_sub)
seurat.integrated.T_cell_sub <- RunPCA(seurat.integrated.T_cell_sub)
ElbowPlot(seurat.integrated.T_cell_sub,ndims = 30)
seurat.integrated.T_cell_sub

seurat.integrated.T_cell_sub <- FindNeighbors(seurat.integrated.T_cell_sub, dims = 1:20)
seurat.integrated.T_cell_sub <- FindClusters(seurat.integrated.T_cell_sub, resolution = c(0.6,0.8,1.0))
seurat.integrated.T_cell_sub <- RunUMAP(seurat.integrated.T_cell_sub, dims = 1:20)
seurat.integrated.T_cell_sub <- RunTSNE(seurat.integrated.T_cell_sub, dims = 1:20)
DimPlot(seurat.integrated.T_cell_sub,reduction = 'umap',group.by='TLS')
DimPlot(seurat.integrated.T_cell_sub,reduction = 'umap',group.by='sample')
save(seurat.integrated.T_cell_sub, file = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/seurat.integrated.B_sub.RData")

B01_Plasma_IGHG1=c(7,10,11,15)
B02_NaiveB_TCL1A=c(1,2,3,9,13,14)
B03_MemoryB_CD27=c(0,6)
B04_GCB_RGS13=c(12)
B05_B_ITGB1=c(4,8)
B06_B_CD3E=c(5)

current.cluster.ids <- c(
  B01_Plasma_IGHG1,
  B02_NaiveB_TCL1A,
  B03_MemoryB_CD27,
  B04_GCB_RGS13,
  B05_B_ITGB1,
  B06_B_CD3E
)

new.cluster.ids <- c(rep("B01_Plasma-IGHG1",length(B01_Plasma_IGHG1)),
                     rep("B02_NaiveB-TCL1A",length(B02_NaiveB_TCL1A)),
                     rep("B03_MemoryB-CD27",length(B03_MemoryB_CD27)),
                     rep("B04_GCB-RGS13",length(B04_GCB_RGS13)),
                     rep("B05_B-ITGB1",length(B05_B_ITGB1)),
                     rep("B06_B-CD3E",length(B06_B_CD3E)))


seurat.integrated.T_cell_sub@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(seurat.integrated.T_cell_sub@meta.data$`integrated_snn_res.1`)), from = current.cluster.ids, to = new.cluster.ids)
table(seurat.integrated.T_cell_sub@meta.data$celltype)
plotCB=as.data.frame(seurat.integrated.T_cell_sub@meta.data)[,"cell"]
DimPlot(seurat.integrated.T_cell_sub, reduction = "umap", group.by = "celltype", pt.size=0.3,cells = plotCB,label=TRUE,label.size = 3)
ggsave(filename = "/NFS2_home/NFS2_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/umap_refined_cluster.pdf", width = 30, height = 20, units = 'cm')
DimPlot(seurat.integrated.T_cell_sub, reduction = "tsne", group.by = "celltype", pt.size=0.3,cells = plotCB,label=TRUE,label.size = 3)
ggsave(filename = "/NFS2_home/NFS2_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/tsne_refined_cluster.pdf", width = 30, height = 20, units = 'cm')
save(seurat.integrated.T_cell_sub,file = "/NFS2_home/NFS2_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/seurat.integrated.B_refined_cluster.Rdata") #Save the Seurat object so it can be loaded directly next time.

# B cell subtype marker gene heatmap
load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/seurat.integrated.B_refined_cluster.Rdata")
Idents(seurat.integrated.T_cell_sub)=seurat.integrated.T_cell_sub@meta.data$celltype
expr=AverageExpression(seurat.integrated.T_cell_sub,assays = "RNA",slot="data")[[1]]
expr <- expr[rowSums(expr)>0,]
expr=as.data.frame(expr)
expr=expr[,c("B01_Plasma-IGHG1","B02_NaiveB-TCL1A","B03_MemoryB-CD27","B04_GCB-RGS13","B05_B-ITGB1","B06_B-CD3E")]
data_top_10=read.table("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/B_subclass_marker_refined_top10.csv",header=T,sep=",")
data_top_10_sort=data_top_10[order(factor(data_top_10$cluster,levels=c(c("B01_Plasma-IGHG1","B02_NaiveB-TCL1A","B03_MemoryB-CD27","B04_GCB-RGS13","B05_B-ITGB1","B06_B-CD3E")))),] 
datamean_plot=expr[data_top_10_sort$gene,]

library(pheatmap)
pheatmap(datamean_plot, 
         show_rownames=T, scale = "row",#show row names (gene names)
         show_colnames=T,
         cluster_cols=F, #do not cluster columns
         cluster_rows=F,
         treeheight_row = 1, #adjust row dendrogram height
         #annotation_col=annotation_col,
         #annotation_row=annotation_row,#add column grouping information into column annotations
         clustering_distance_rows = "correlation", #optimize row clustering distance
         color = colorRampPalette(c("#FF00FF", "#000000", "#FFFF00"))(10000),
         #annotation_colors=ann_colors, #assign colors to annotation blocks
         #cutree_rows=3, #split all rows into three groups based on the clustering tree
         gaps_col=c(1:6), #set column gap positions
         #gaps_row=c(7,13),
         cellheight=4,cellwidth=8, #size of each heatmap cell
         fontsize_col=4, #font size for column names
         fontsize_row = 4, #font size for row names
         fontsize_annotation=0.1,
         angle_col=315,
         border=FALSE,
         main = "Top 10 of B_cell", #main title
         filename="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/marker_avg_refined_top10.png")

# B cell monocle analysis
library(monocle)

load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/seurat.integrated.B_refined_cluster.Rdata")
DefaultAssay(seurat.integrated.T_cell_sub)="RNA"
expr_matrix <- as(as.matrix(seurat.integrated.T_cell_sub@assays$RNA@counts), 'sparseMatrix')
p_data <- seurat.integrated.T_cell_sub@meta.data
f_data <- data.frame(gene_short_name = row.names(seurat.integrated.T_cell_sub),row.names = row.names(seurat.integrated.T_cell_sub))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#Convert p_data and f_data from data.frame to AnnotatedDataFrame objects.
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.5) #This step adds a num_cells_expressed column to fData(cds).
print(head(fData(cds)))#There are 13,714 genes at this stage.
#expressed_genes <- row.names(subset(fData(cds), num_cells_expressed > 10))
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed > nrow(pd) * 0.01))
test_cds=cds
diff <- differentialGeneTest(test_cds[expressed_genes,],fullModelFormulaStr="~celltype",cores=10) 
deg <- subset(diff, qval < 1e-5) #Select 2,829 genes.
deg <- deg[order(deg$qval,decreasing=F),]
ordergene <- row.names(deg)[order(deg$qval)][1:500]
test_cds <- setOrderingFilter(cds, ordergene)  
#This is a critical step: after selecting genes, use setOrderingFilter to store them in the cds object.
#After setOrderingFilter, selected genes are stored in cds@featureData@data[["use_for_ordering"]], which can be inspected with table(...).
#pdf("cluster_visualization_B/train.ordergenes.pdf")
plot_ordering_genes(test_cds)
#dev.off()
#In this plot, black points are trajectory-building genes, gray points are background genes, and the red curve is the trend from step 2.
test_cds <- reduceDimension(test_cds, max_components = 2,
    method = 'DDRTree')

#test_cds <- orderCells(test_cds)
test_cds <- orderCells(test_cds,root_state = 5)
save(test_cds,file="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/test_cds.Rdata")
save(ordergene,file="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/ordergene.Rdata")

# plot trajectory
p=plot_cell_trajectory(test_cds,color_by="Pseudotime",cell_size=0.5, show_tree=FALSE,show_backbone=FALSE,show_branch_points=FALSE) 
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/trajecory_pseudotime.pdf", plot = p, device = 'pdf', width = 15, height = 15, units = 'cm')
p=plot_cell_trajectory(test_cds,color_by = "celltype",cell_size=0.5, show_tree=FALSE,show_backbone=FALSE,show_branch_points=FALSE)+scale_color_manual(values=c("#64c6d7","#4579ba","#ec9730","#e74b4e","#f2c844","#c667a4"))
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/trajecory_celltype.pdf", plot = p, device = 'pdf', width = 15, height = 15, units = 'cm')
p=plot_cell_trajectory(test_cds,color_by = "celltype",cell_size=0.5, show_tree=FALSE,show_backbone=FALSE,show_branch_points=FALSE)+scale_color_manual(values=c("#64c6d7","#4579ba","#ec9730","#e74b4e","#f2c844","#c667a4"))
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/trajecory_celltype.pdf", plot = p, device = 'pdf', width = 15, height = 15, units = 'cm')
p=plot_cell_trajectory(test_cds,color_by = "State")
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/trajecory_state.png", plot = p, device = 'png', width = 30, height = 30, units = 'cm')
p=plot_cell_trajectory(test_cds,color_by = "celltype",cell_size=0.5, show_tree=FALSE,show_backbone=FALSE,show_branch_points=FALSE)+facet_wrap(~celltype,nrow=1)
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/trajecory_facet_celltype.png", plot = p, device = 'png', width = 40, height = 10, units = 'cm')
p=plot_cell_trajectory(test_cds,color_by = "celltype",cell_size=0.5, show_tree=FALSE,show_backbone=FALSE,show_branch_points=FALSE)+scale_color_manual(values=c("#64c6d7","#4579ba","#ec9730","#e74b4e","#f2c844","#c667a4"))+facet_wrap(~TLS,nrow=2)
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/trajecory_facet_TLS.pdf", plot = p, device = 'pdf', width = 15, height = 30, units = 'cm')
p=plot_cell_trajectory(test_cds,color_by = "celltype")+facet_wrap(~sample,nrow=1)
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/trajecory_facet_sample.png", plot = p, device = 'png', width = 40, height = 10, units = 'cm')

### Pseudotime density plot
library(ggpubr)
df <- pData(test_cds) 
library(ggridges)
p=ggplot(df, aes(x=Pseudotime,y=TLS,fill=TLS))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,7.5),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )+scale_color_manual(values=c("#2a9594","#c667a4"))

ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/Pseudotime_density_TLS.pdf", plot = p, device = 'png', width = 40, height = 10, units = 'cm')

p=ggplot(df, aes(x=Pseudotime,y=celltype,fill=celltype))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,7.5),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/Pseudotime_density.pdf", plot = p, device = 'png', width = 40, height = 10, units = 'cm')


# Cluster tree
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(6,"Dark2")
colorl <- rep(color,each=10)
library("ggtree")
load(paste0("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/seurat.integrated.B_refined_cluster.Rdata"))
Idents(seurat.integrated.T_cell_sub)=seurat.integrated.T_cell_sub@meta.data$celltype
color <- brewer.pal(6,"Dark2")
seurat <- BuildClusterTree(
  seurat.integrated.T_cell_sub,
  dims = 1:15,
  reorder = FALSE,
  reorder.numeric = FALSE
)

tree <- seurat@tools$BuildClusterTree
pdf(file='/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/cluster_tree.pdf', height = 4, width = 8)
ggtree::ggtree(tree, aes(x, y)) +
  scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(color = color, shape = 16, size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,3,0,0), 'cm'))
dev.off()


TLS=seurat.integrated.T_cell_sub[,which(seurat.integrated.T_cell_sub@meta.data$TLS=="With_TLS")]
NoTLS=seurat.integrated.T_cell_sub[,which(seurat.integrated.T_cell_sub@meta.data$TLS=="No_TLS")]
seurat <- BuildClusterTree(
  TLS,
  dims = 1:15,
  reorder = FALSE,
  reorder.numeric = FALSE
)

tree <- seurat@tools$BuildClusterTree

pdf(file='/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/cluster_tree_TLS.pdf', height = 4, width = 8)

ggtree::ggtree(tree, aes(x, y)) +
  scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(color = color, shape = 16, size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,3,0,0), 'cm'))

dev.off()
seurat <- BuildClusterTree(
  NoTLS,
  dims = 1:15,
  reorder = FALSE,
  reorder.numeric = FALSE
)

tree <- seurat@tools$BuildClusterTree

pdf(file='/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/cluster_tree_NoTLS.pdf', height = 4, width = 8)
ggtree::ggtree(tree, aes(x, y)) +
  scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(color = color, shape = 16, size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,3,0,0), 'cm'))

dev.off()

# Cluster distance
library(reshape2)
load(paste0("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/seurat.integrated.B_refined_cluster.Rdata"))
Idents(seurat.integrated.T_cell_sub)=seurat.integrated.T_cell_sub@meta.data$celltype
pca_embeddings=seurat.integrated.T_cell_sub@reductions$pca@cell.embeddings
E_Distance_matrix <- data.frame(celltypes1=character(),celltypes2=character(),Distance=double())

for (i in sort(unique(seurat.integrated.T_cell_sub@meta.data$celltype))){
    for (j in sort(unique(seurat.integrated.T_cell_sub@meta.data$celltype))){
        tmp_name <- paste0(i,'_',j)
        type1=seurat.integrated.T_cell_sub[,which(seurat.integrated.T_cell_sub@meta.data$celltype==i)]
        type2=seurat.integrated.T_cell_sub[,which(seurat.integrated.T_cell_sub@meta.data$celltype==j)]
        type1_num=length(type1@meta.data$celltype)
        type2_num=length(type2@meta.data$celltype)
        min_num=min(type1_num,type2_num)
        select_num=50
        set.seed(40)
        type1_select = type1[, sample(1:type1_num,select_num)]
        set.seed(40)
        type2_select = type2[, sample(1:type2_num,select_num)]
        type1_select_pca=type1_select$pca@cell.embeddings
        type2_select_pca=type2_select$pca@cell.embeddings
        E_Distance=median(rowSums((type1_select_pca-type2_select_pca)**2)^0.5)
        E_Distance_matrix[tmp_name,"celltypes1"]=i
        E_Distance_matrix[tmp_name,"celltypes2"]=j
        E_Distance_matrix[tmp_name,"Distance"]=E_Distance
    }
}

write.table(E_Distance_matrix,file="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/E_Distance_matrix.txt",quote = F,sep="\t",row.names = F)

data_dis=read.table("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/E_Distance_matrix.txt",header = T,sep="\t")
data_dis_dcast=dcast(data_dis,celltypes1~celltypes2)
rownames(data_dis_dcast)=data_dis_dcast$celltypes1
data_dis_dcast$celltypes1=NULL
data_dis_dcast=as.matrix(data_dis_dcast)

library(pheatmap)
pheatmap(data_dis_dcast, 
         show_rownames=T, scale = "none",#show row names (gene names)
         show_colnames=T,
         cluster_cols=F, #do not cluster columns
         cluster_rows=F,
         treeheight_row = 1, #adjust row dendrogram height
         #annotation_col=annotation_col,
         #annotation_row=annotation_row,#add column grouping information into column annotations
         clustering_distance_rows = "correlation", #optimize row clustering distance
         color = colorRampPalette(c(color <- brewer.pal(6,"Accent")[4:6]))(10000),
         #annotation_colors=ann_colors, #assign colors to annotation blocks
         #cutree_rows=3, #split all rows into three groups based on the clustering tree
         #gaps_col=c(1:22), #set column gap positions
         #gaps_row=c(7,13),
         cellheight=30,cellwidth=30, #size of each heatmap cell
         fontsize_col=10, #font size for column names
         fontsize_row = 10, #font size for row names
         fontsize_annotation=0.1,
         angle_col=315,
         border=FALSE,
         filename="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/E_distance.png")

# CytoTRACE
library("CytoTRACE")
suppressMessages(library(Seurat))
suppressMessages(library(DoubletFinder))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(clusterProfiler))
library(msigdbr)
library(velocyto.R)
library(DoubletFinder)
load("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/seurat.integrated.B_refined_cluster.Rdata")

count_matrix=as.matrix(seurat.integrated.T_cell_sub@assays$RNA@counts)
metadata=as.character(seurat.integrated.T_cell_sub@meta.data$celltype)
names(metadata) <- rownames(seurat.integrated.T_cell_sub@meta.data)

results <- CytoTRACE(mat = count_matrix,ncores=10)
embedding=seurat.integrated.T_cell_sub@reductions$umap@cell.embeddings
plotCytoGenes(results, numOfGenes = 50,outputDir = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/CytoTRACE")
plotCytoTRACE(results, phenotype = metadata,emb=embedding ,outputDir = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_B/CytoTRACE")
# CytoTRACE score and pseudotime
p=ggplot(B_meta, aes(Pseudotime, CytoTRACE)) +
  geom_smooth ()+
ylab("CytoTRACE score")+
labs(title="CytoTRACE score")+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    plot.title = element_text(size=20,hjust=0.5,color="black"),
    axis.text.x = element_text(hjust=1,vjust=0.6),
    #legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    #legend.direction = "horizontal",
    #legend.position = c(0.5,0.9),
    legend.background = element_blank(),
    #panel.background = element_rect(fill = "transparent",colour = "black"),
    plot.background = element_blank()
  )
# T/NK cell recluster and annotation

seurat.integrated.bigCluster=readRDS("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_bigcluster/seurat.integrated.bigCluster.rds")

scRNA_T_cell <- seurat.integrated.bigCluster[,which(seurat.integrated.bigCluster@meta.data$celltype=="T_cell" | seurat.integrated.bigCluster@meta.data$celltype=="NK_cell")]
scRNA_T_cell_filter <- subset(scRNA_T_cell, subset = nFeature_RNA > 500  & percent.mt < 10 & nCount_RNA>200 & log10GenesPerUMI > 0.80)
scRNA_T_cell_sub=CreateSeuratObject(counts = scRNA_T_cell_filter@assays$RNA@counts,
                       meta.data = scRNA_T_cell_filter@meta.data)

obj.list <- SplitObject(scRNA_T_cell_sub, split.by = 'sample')
for(i in 1:length(obj.list)){
    obj.list[[i]] <- NormalizeData(object=obj.list[[i]])
    obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}
features <- SelectIntegrationFeatures(object.list=obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                 anchor.features = features)
seurat.integrated.T_cell_sub <- IntegrateData(anchorset = anchors)
seurat.integrated.T_cell_sub <- ScaleData(seurat.integrated.T_cell_sub)
seurat.integrated.T_cell_sub <- RunPCA(seurat.integrated.T_cell_sub)
ElbowPlot(seurat.integrated.T_cell_sub,ndims = 30)
seurat.integrated.T_cell_sub <- FindNeighbors(seurat.integrated.T_cell_sub, dims = 1:20)
seurat.integrated.T_cell_sub <- FindClusters(seurat.integrated.T_cell_sub, resolution = c(0.4,0.6,0.8,1))
seurat.integrated.T_cell_sub <- RunUMAP(seurat.integrated.T_cell_sub, dims = 1:20)
seurat.integrated.T_cell_sub <- RunTSNE(seurat.integrated.T_cell_sub, dims = 1:20)
save(seurat.integrated.T_cell_sub, file = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_T/seurat.integrated.T_sub.RData")

L01_NK_TRDC=c(8)
L02_NK_KIR2DL3=c(13)
L03_NK_FGFBP2=c(6)
L04_NK_GNLY=c(19)
L05_ILC_KIT=c(24)
T01_CD4_CD40LG=c(10)
T02_CD4_TCF7=c(1)  
T03_CD4_CCR6=c(15)  
T04_CD4_LTB=c(9)  
T05_CD4_FOXP3=c(4,21,23,25)
T06_CD4_CXCL13=c(18)
T07_CD8_PDCD1=c(3)
T08_CD8_ZNF683=c(0)
T09_CD8_KLF2=c(7)
T10_CD8_GZMA=c(5)
T11_CD8_APOE=c(16)
T12_CD8_ISG=c(11)
T13_CD8_STMN1=c(17)
T14_CD8_MKI67=c(12)
T15_T_TTN=c(2)
T16_gdT_TRDV2=c(20)
T17_MAIT_SLC4A10=c(14)




current.cluster.ids <- c(
L01_NK_TRDC,
L02_NK_KIR2DL3,
L03_NK_FGFBP2,
L04_NK_GNLY,
L05_ILC_KIT,
T01_CD4_CD40LG,
T02_CD4_TCF7,
T03_CD4_CCR6,
T04_CD4_LTB,
T05_CD4_FOXP3,
T06_CD4_CXCL13,
T07_CD8_PDCD1,
T08_CD8_ZNF683,
T09_CD8_KLF2,
T10_CD8_GZMA,
T11_CD8_APOE,
T12_CD8_ISG,
T13_CD8_STMN1,
T14_CD8_MKI67,
T15_T_TTN,
T16_gdT_TRDV2,
T17_MAIT_SLC4A10
)

new.cluster.ids <- c(rep("L01_NK-TRDC",length(L01_NK_TRDC)),
                     rep("L02_NK-KIR2DL3",length(L02_NK_KIR2DL3)),
                     rep("L03_NK-FGFBP2",length(L03_NK_FGFBP2)),
                     rep("L04_NK-GNLY",length(L04_NK_GNLY)),
                     rep("L05_ILC-KIT",length(L05_ILC_KIT)),
                     rep("T01_CD4-CD40LG",length(T01_CD4_CD40LG)),
                     rep("T02_CD4-TCF7",length(T02_CD4_TCF7)),
                     rep("T03_CD4-CCR6",length(T03_CD4_CCR6)),
                     rep("T04_CD4-LTB",length(T04_CD4_LTB)),
                     rep("T05_CD4-FOXP3",length(T05_CD4_FOXP3)),
                     rep("T06_CD4-CXCL13",length(T06_CD4_CXCL13)),
                     rep("T07_CD8-PDCD1",length(T07_CD8_PDCD1)),
                     rep("T08_CD8-ZNF683",length(T08_CD8_ZNF683)),
                     rep("T09_CD8-KLF2",length(T09_CD8_KLF2)),
                     rep("T10_CD8-GZMA",length(T10_CD8_GZMA)),
                     rep("T11_CD8-APOE",length(T11_CD8_APOE)),
                     rep("T12_CD8-ISG",length(T12_CD8_ISG)),
                     rep("T13_CD8-STMN1",length(T13_CD8_STMN1)),
                     rep("T14_CD8-MKI67",length(T14_CD8_MKI67)),
                     rep("T15_T-TTN",length(T15_T_TTN)),
                     rep("T16_gdT-TRDV2",length(T16_gdT_TRDV2)),
                     rep("T17_MAIT-SLC4A10",length(T17_MAIT_SLC4A10))
                     
)

seurat.integrated.T_cell_sub@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(seurat.integrated.T_cell_sub@meta.data$`integrated_snn_res.1`)), from = current.cluster.ids, to = new.cluster.ids)
table(seurat.integrated.T_cell_sub@meta.data$celltype)
plotCB=as.data.frame(seurat.integrated.T_cell_sub@meta.data)[,"cell"]
DimPlot(seurat.integrated.T_cell_sub, reduction = "umap", group.by = "celltype", pt.size=0.3,cells = plotCB,label=TRUE,label.size = 3)
ggsave(filename = "results/cluster_visualization_T/umap_refined_cluster.pdf", width = 40, height = 24, units = 'cm')
DimPlot(seurat.integrated.T_cell_sub, reduction = "tsne", group.by = "celltype", pt.size=0.3,cells = plotCB,label=TRUE)
ggsave(filename = "results/cluster_visualization_T/tsne_refined_cluster.pdf", width = 40, height = 24, units = 'cm')

save(seurat.integrated.T_cell_sub,file = "results/cluster_visualization_T/seurat.integrated.T_refined_cluster.Rdata")


# Myeloid recluster
seurat.integrated.bigCluster=readRDS("/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_bigcluster/seurat.integrated.bigCluster.rds")
scRNA_T_cell <- seurat.integrated.bigCluster[,which(seurat.integrated.bigCluster@meta.data$celltype=="Myeloid")]
scRNA_T_cell_filter <- subset(scRNA_T_cell, subset = nFeature_RNA > 500  & percent.mt < 10 & nCount_RNA>200 & log10GenesPerUMI > 0.80)

scRNA_T_cell_sub=CreateSeuratObject(counts = scRNA_T_cell_filter@assays$RNA@counts,
                       meta.data = scRNA_T_cell_filter@meta.data)

obj.list <- SplitObject(scRNA_T_cell_sub, split.by = 'sample')
for(i in 1:length(obj.list)){
    obj.list[[i]] <- NormalizeData(object=obj.list[[i]])
    obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}
features <- SelectIntegrationFeatures(object.list=obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                 anchor.features = features)
seurat.integrated.T_cell_sub <- IntegrateData(anchorset = anchors)

seurat.integrated.T_cell_sub <- ScaleData(seurat.integrated.T_cell_sub)
seurat.integrated.T_cell_sub <- RunPCA(seurat.integrated.T_cell_sub)
ElbowPlot(seurat.integrated.T_cell_sub,ndims = 30)

seurat.integrated.T_cell_sub <- FindNeighbors(seurat.integrated.T_cell_sub, dims = 1:20)
seurat.integrated.T_cell_sub <- FindClusters(seurat.integrated.T_cell_sub, resolution = c(0.4,0.6,0.8,1.0))
seurat.integrated.T_cell_sub <- RunUMAP(seurat.integrated.T_cell_sub, dims = 1:20)
seurat.integrated.T_cell_sub <- RunTSNE(seurat.integrated.T_cell_sub, dims = 1:20)
save(seurat.integrated.T_cell_sub, file = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_myeloid/seurat.integrated.Myeloid_sub.RData")

M01_mDC_LAMP3=c(22)
M02_cDC1_BATF3=c(20)
M03_cDC2_CD1C=c(3)
M04_Mono_S100A8=c(6)
M05_Mono_FCGR3A=c(17)
M06_Mono_FN1=c(4,10)
M07_Mono_CD3E=c(7,15,18)
M08_Macro_MT1M=c(16)
M09_Macro_STAB1=c(9)
M10_Macro_CCL3L1=c(1)
M11_Macro_IFI27=c(5)
M12_Macro_OLFML3=c(11)
M13_Macro_CXCL10=c(8)
M14_Macro_FOLR2=c(2,19)
M15_Macro_SPP1=c(0)
M16_Macro_CCL18=c(13)
M17_Macro_STMN1=c(12,14)


current.cluster.ids <- c(
M01_mDC_LAMP3,
M02_cDC1_BATF3,
M03_cDC2_CD1C,
M04_Mono_S100A8,
M05_Mono_FCGR3A,
M06_Mono_FN1,
M07_Mono_CD3E,
M08_Macro_MT1M,
M09_Macro_STAB1,
M10_Macro_CCL3L1,
M11_Macro_IFI27,
M12_Macro_OLFML3,
M13_Macro_CXCL10,
M14_Macro_FOLR2,
M15_Macro_SPP1,
M16_Macro_CCL18,
M17_Macro_STMN1
)

new.cluster.ids <- c(rep("M01_mDC-LAMP3",length(M01_mDC_LAMP3)),
                     rep("M02_cDC1-BATF3",length(M02_cDC1_BATF3)),
                     rep("M03_cDC2-CD1C",length(M03_cDC2_CD1C)),
                     rep("M04_Mono-S100A8",length(M04_Mono_S100A8)),
                     rep("M05_Mono-FCGR3A",length(M05_Mono_FCGR3A)),
                     rep("M06_Mono-FN1",length(M06_Mono_FN1)),
                     rep("M07_Mono-CD3E",length(M07_Mono_CD3E)),
                     rep("M08_Macro-MT1M",length(M08_Macro_MT1M)),
                     rep("M09_Macro-STAB1",length(M09_Macro_STAB1)),
                     rep("M10_Macro-CCL3L1",length(M10_Macro_CCL3L1)),
                     rep("M11_Macro-IFI27",length(M11_Macro_IFI27)),
                     rep("M12_Macro-OLFML3",length(M12_Macro_OLFML3)),
                     rep("M13_Macro-CXCL10",length(M13_Macro_CXCL10)),
                     rep("M14_Macro-FOLR2",length(M14_Macro_FOLR2)),
                     rep("M15_Macro-SPP1",length(M15_Macro_SPP1)),
                     rep("M16_Macro-CCL18",length(M16_Macro_CCL18)),
                     rep("M17_Macro-STMN1",length(M17_Macro_STMN1))
)

seurat.integrated.T_cell_sub@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(seurat.integrated.T_cell_sub@meta.data$`integrated_snn_res.1`)), from = current.cluster.ids, to = new.cluster.ids)
table(seurat.integrated.T_cell_sub@meta.data$celltype)
plotCB=as.data.frame(seurat.integrated.T_cell_sub@meta.data)[,"cell"]
DimPlot(seurat.integrated.T_cell_sub, reduction = "umap", group.by = "celltype", pt.size=0.5,cells = plotCB,label=TRUE)
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_myeloid/umap_refined_cluster.pdf", width = 15, height = 12, units = 'cm')
DimPlot(seurat.integrated.T_cell_sub, reduction = "tsne", group.by = "celltype", pt.size=0.5,cells = plotCB,label=TRUE)
ggsave(filename = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_myeloid/tsne_refined_cluster.pdf", width = 15, height = 12, units = 'cm')
save(seurat.integrated.T_cell_sub,file = "/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_myeloid/seurat.integrated.myeloid_refined_cluster.Rdata") 

DimPlot(seurat.integrated.T_cell_sub,reduction = 'umap',group.by='TLS')
ggsave(filename="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_myeloid/myeloid_subcluster_groupbyTLS.png",width=15,height=15,units='cm')
DimPlot(seurat.integrated.T_cell_sub,reduction = 'umap',group.by='sample')
ggsave(filename="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_myeloid/myeloid_subcluster_groupbySample.png",width=15,height=15,units='cm')


Idents(seurat.integrated.T_cell_sub)=seurat.integrated.T_cell_sub@meta.data$celltype
DimPlot(seurat.integrated.T_cell_sub, reduction = "umap",label=TRUE,split.by = "TLS",label.size = 3)
ggsave(filename="/NFS_home/NFS_home_2/zhouchi/ICC/scRNA/results/cluster_visualization_myeloid/myeloid_subcluster_groupbyTLS_facet.png",width=30,height=15,units='cm')





