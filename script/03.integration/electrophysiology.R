#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("reshape2"))

# * read data
patch_data <- readRDS("data/meta/patch_data.rds")
sample.dict <- readRDS("data/meta/sample.dict.rds")
patch_logcpm_filter <- readRDS("data/meta/patch_logcpm_filter.rds")
identical(colnames(patch_logcpm_filter), as.character(sample.dict$CellID))
rownames(sample.dict) <- sample.dict$CellID
ion_data <- read.table("data/resource/IonChannelGene.txt", sep="\t", header=TRUE)
ion_330genes <- str_to_title(ion_data$Approved.symbol)

# * spearman correlation test
patch_logcpm_ion330genes <- patch_logcpm_filter[which(rownames(patch_logcpm_filter) %in% ion_330genes),]
patch_logcpm_ion330genes <- patch_logcpm_ion330genes[rowSums(patch_logcpm_ion330genes)>0, ]
electro_continuous <- c("Capacitance", "Input.resistance", "Time.constant", "Resting.membrane.potential", "AP.Frequency", 
	"Ramp.AP.number", "Rheobase", "Max.AP.number", "AP.Amplitude", "AP.half.width", 
	"AP.threshold", "AHP")

patch_iongene_test <- list()
for (i in electro_continuous) {
	sample.dict %>% select(i) -> sample.dict_select
	electro_select <- as.numeric(sample.dict_select[,1])
	names(electro_select) <- rownames(sample.dict_select)
	electro_select <- electro_select[!is.na(electro_select)]
	patch_logcpm_select <- patch_logcpm_ion330genes[,match(names(electro_select), colnames(patch_logcpm_ion330genes))]
	if (identical(names(electro_select), colnames(patch_logcpm_select))){
		spearman_test <- data.frame(rho=rep(NA, nrow(patch_logcpm_select)), pval=rep(NA, nrow(patch_logcpm_select)))
		rownames(spearman_test) <- rownames(patch_logcpm_select)
		for (j in 1:nrow(patch_logcpm_select)){
			spearman_test[j, "rho"] <- cor.test(as.numeric(electro_select), as.numeric(patch_logcpm_select[j,]), method="spearman")$estimate
			spearman_test[j, "pval"] <- cor.test(as.numeric(electro_select), as.numeric(patch_logcpm_select[j,]), method="spearman")$p.value
		}
		spearman_test$padj <- p.adjust(spearman_test$pval, method="BH")
		patch_iongene_test[[i]] <- spearman_test 
	}
}

# * two-sided unpaired T test
electro_discrete <- c("Rebound", "Burst.AP", "Spont.action.potentials")
for (i in electro_discrete) {
	sample.dict %>% select(i) -> sample.dict_select
	electro_select <- as.factor(sample.dict_select[,1])
	names(electro_select) <- rownames(sample.dict_select)
	electro_select <- electro_select[which(electro_select!="NaN")]
	patch_logcpm_select <- patch_logcpm_ion330genes[,match(names(electro_select), colnames(patch_logcpm_ion330genes))]
	if (identical(names(electro_select), colnames(patch_logcpm_select))){
		ttest_test <- data.frame(log2FC=rep(NA, nrow(patch_logcpm_select)), pval=rep(NA, nrow(patch_logcpm_select)))
		rownames(ttest_test) <- rownames(patch_logcpm_select)
		for (j in 1:nrow(patch_logcpm_select)){
			ttest_test[j, "log2FC"] <- log2(mean(as.numeric(patch_logcpm_select[j,which(electro_select==1)]))/(mean(as.numeric(patch_logcpm_select[j,which(electro_select==0)]))+1e-10))
			ttest_test[j, "pval"] <- wilcox.test(as.numeric(patch_logcpm_select[j,which(electro_select==1)]), as.numeric(patch_logcpm_select[j,which(electro_select==0)]))$p.value
		}
		ttest_test$padj <- p.adjust(ttest_test$pval, method="BH")
		patch_iongene_test[[i]] <- ttest_test 
	}
}

# * output
write.csv(patch_iongene_test, file="data/meta/electro_Ion312Gene_test.csv")

