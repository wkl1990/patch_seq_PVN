#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("dendsort"))

# * read data
PVH_count <- read.csv("data/resource/GSE148568/GSE148568_compiled_data.csv", row.names=1)
PVH_meta <- read.csv("data/resource/GSE148568/GSE148568_cell_metadata_after_qc.csv", row.names=1)
gene_count_mtx <- readRDS("rds/gene_count_mtx.rds")
gene_cpm <- readRDS("rds/gene_cpm.rds")
PVH_count_filter <- PVH_count[,match(PVH_meta$sample_id, colnames(PVH_count))]
patch_cpm_filter <- gene_cpm[,which(colnames(gene_cpm) %in% sample_information_filter$CellID)]

# marker gene list
markergene_list <- c("Sim1", "Slc17a6", "Slc32a1", "Gad2", "Oxt", "Trac", "Cela1", "Pdyn", "Avp", "C1ql2", 
	"Ebf3", "Igfbp2", "Mfge8", "Onecut2", "Ret", "Six3", "Coch", "Dlx1", "Dlx2", "Dlx6", 
	"Npy1r", "Aldh1a1", "Col12a1", "Kcnv1", "Sfrp2", "Galr1", "Sst", "Kcnh8", "Gpr101", "Klf4",
	"Fign", "Chodl", "Trh", "Pnoc", "Ebf1", "Crh", "Scgn", "Bcl11a", "Gbx2", "Grm2", 
	"Lhx9", "Kitl", "Inhbb", "Sox5", "Brs3", "Cdkn1a", "Fam150b", "Calb2", "Ntng1", "Reln", 
	"Penk")
gene_list <- c("GAD1", "GAD2", "SLC17A7", "SLC17A6", "SLC17A8", "SLC32A1", "GLP1R", "GLP2R", "MC4R", "CRFR1", 
	"CRFR2", "LEPR", "SIM1", "Npy1r", "Crh", "Reln", "Ntng1", "Pdyn", "Penk", "TRH", "Oxt", "Avp", "SST", "ADCYAP", "BDNF", "SYT1", "SYT2")
gene_list <- str_to_title(gene_list)
markergene_list <- unique(c(markergene_list,gene_list)) 
overlap_genes <- intersect(rownames(gene_cpm_filter),rownames(PVH_count_filter))
overlap_markergenes <- intersect(markergene_list,overlap_genes)

# * clustering and cell annotation
PVH_cpm <- cpm(PVH_count_filter)
PVH_logcpm <- log2(PVH_cpm + 1)
PVH_logcpm_markergenes <- PVH_logcpm[match(overlap_markergenes, rownames(PVH_logcpm)), ]
PVH_zscorelogcpmsample_markergenes <- scale(PVH_logcpm_markergenes)
patch_logcpm <- log2(patch_cpm_filter + 1)
patch_logcpm_markergenes <- patch_logcpm[match(overlap_markergenes, rownames(patch_logcpm)), ]
patch_zscorelogcpmsample_markergenes <- scale(patch_logcpm_markergenes)
if (identical(rownames(patch_zscorelogcpmsample_markergenes), rownames(PVH_zscorelogcpmsample_markergenes))){
	combined_zscorelogcpmsample_markergenes <- cbind(PVH_zscorelogcpmsample_markergenes, patch_zscorelogcpmsample_markergenes)
} 
combined_zscorelogcpmsample_markergenes_combat <- ComBat(dat=combined_zscorelogcpmsample_markergenes, batch=c(rep(1,ncol(PVH_zscorelogcpmsample_markergenes)), rep(2,ncol(patch_zscorelogcpmsample_markergenes))), mod=NULL, par.prior=TRUE, ref.batch=1)
combined_zscorelogcpmsample_dist <- dist(t(combined_zscorelogcpmsample_markergenes_combat), method="euclidean")
combined_zscorelogcpmsample_hclust <- hclust(combined_zscorelogcpmsample_dist, method="ward.D2")
patch_celltype_clust <- cutreeStatic(combined_zscorelogcpmsample_hclust, cutHeight=60, minSize=10)
names(patch_celltype_clust) <- combined_zscorelogcpmsample_hclust$labels
Oxt_cells <- patch_celltype_clust[patch_celltype_clust==2&grepl("RNA", names(patch_celltype_clust))]
Other_cells <- patch_celltype_clust[patch_celltype_clust!=2&grepl("RNA", names(patch_celltype_clust))]

# * save cell annotation
saveRDS(Oxt_cells, file="data/meta/Oxt_cells.rds")
saveRDS(Other_cells, file="data/meta/Other_cells.rds")


