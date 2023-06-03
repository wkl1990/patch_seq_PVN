#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("edgeR"))

# * read data
gene_count_mtx <- readRDS(file="data/meta/gene_count_mtx.rds")
gene_cpm <- readRDS(file="data/meta/gene_cpm.rds")
sample_information <- readRDS(file="data/meta/sample_information.rds")


# * differential gene expression 
sample_information_filter <- sample_information %>% filter(QC=="Pass")
gene_count_filter <- gene_count_mtx[,which(colnames(gene_count_mtx) %in% sample_information_filter$CellID)]
design_filter <- data.frame(
    SampleID=sample_information_filter$Sample.Number[match(colnames(gene_count_filter), sample_information_filter$CellID)],
    condition=sample_information_filter$condition[match(colnames(gene_count_filter), sample_information_filter$CellID)],
    QC=sample_information_filter$QC[match(colnames(gene_count_filter), sample_information_filter$CellID)], 
    stringsAsFactors=TRUE)
gene_dataset_filter <- DESeqDataSetFromMatrix(gene_count_filter, DataFrame(design_filter), ~condition)
# keep only rows that have at least 1 reads total
keep <- rowSums(counts(gene_dataset_filter)) >= 1
gene_dataset_filter <- gene_dataset_filter[keep,]
# set control group
gene_dataset_filter$condition <- relevel(gene_dataset_filter$condition, ref="Fed")
gene_dataset_filter <- DESeq(gene_dataset_filter)
resLFC <- lfcShrink(gene_dataset_filter, coef="condition_Fasted_vs_Fed", type="ashr")
resLFC.Ordered <- resLFC[order(resLFC$pvalue),]
resLFC.Ordered %>% as.data.frame %>% filter(!is.na(padj)) -> resLFC.Ordered.nonNA
write.csv(as.data.frame(resLFC.Ordered.nonNA), file="data/meta/Fasted_vs_Fed.LFC.nonNA.csv")

