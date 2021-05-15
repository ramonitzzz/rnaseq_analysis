#setwd("~/OneDrive/Bioinfo/Proyectos/Tuberculosis/RNA-seq/Run5/Fernando")
#setwd("~/Desktop/20_RNAseq/Adiana/Bovis_only")
setwd('/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes')



#install.packages("readxl")

library("readxl")
library("xlsx")
library(goseq)
library(randomcoloR)
library(ggplot2)
library(Rsubread)
library(edgeR)
library(statmod)
library(calibrate)
library(gplots)
library(PCAtools)

dir.create("heatmapfiles")
dir.create("Pictures")
dir.create("Pictures/VolcanoA")
dir.create("Pictures/VolcanoB")
dir.create("Tables")
dir.create("Tables/FDR")
dir.create("Tables/GOseq/")
dir.create("Tables/GOseq/bp/")
dir.create("Tables/GOseq/mf/")
dir.create("Tables/GOseq/cc/")
dir.create("Pictures/GEA/")
dir.create("Pictures/GEA/bp/")
dir.create("Pictures/GEA/mf/")
dir.create("Pictures/GEA/cc/")
dir.create("Pictures/GEA/bp/")
dir.create("Pictures/GEA/bp/piechart/")
dir.create("Pictures/GEA/mf/piechart/")
dir.create("Pictures/GEA/cc/piechart/")

# Load data and set conditions
samplesinfo <- c("IB-p0;pH-6-3_PZA-0", "IB-p50;pH-6-3_PZA-50",
                 "IC-p0;pH-6-3_PZA-0", "IC-p50;pH-6-3_PZA-50",
                 "IIB-p0;pH-6-3_PZA-0", "IIB-p50;pH-6-3_PZA-50",
                 "IIC-p0;pH-6-3_PZA-0", "IIC-p50;pH-6-3_PZA-50",
                 "IIIB-p0;pH-6-3_PZA-0", "IIIB-p50;pH-6-3_PZA-50",
                 "IIIC-p0;pH-6-3_PZA-0", "IIIC-p50;pH-6-3_PZA-50")

# Try to group by treat + sample

Treat2 <- factor(c("Bwt-p0;pH-6-3_PZA-0", "Bwt-p50;pH-6-3_PZA-50",
                   "Bpnca-p0;pH-6-3_PZA-0", "Bpnca-p50;pH-6-3_PZA-50",
                   "Bwt-p0;pH-6-3_PZA-0", "Bwt-p50;pH-6-3_PZA-50",
                   "Bpnca-p0;pH-6-3_PZA-0", "Bpnca-p50;pH-6-3_PZA-50",
                   "Bwt-p0;pH-6-3_PZA-0", "Bwt-p50;pH-6-3_PZA-50",
                   "Bpnca-p0;pH-6-3_PZA-0", "Bpnca-p50;pH-6-3_PZA-50"))


samplesgroup <-c("Bwt0", "Bwt50", "Bpnca0", "Bpnca50",
                 "Bwt0", "Bwt50", "Bpnca0", "Bpnca50",
                 "Bwt0", "Bwt50", "Bpnca0", "Bpnca50")


# Counts
counts_salmon = read.table("counts_bovis.tsv", sep="\t", header = TRUE, row.names = 1, stringsAsFactors = F, check.names = T)

#fc.dge <- DGEList(counts=fc$counts, genes = fc$annotation, samples = samplesinfo, group = samplesgroup)
fc.dge <- DGEList(counts = counts_salmon, genes = row.names(counts_salmon), 
                  samples = colnames(counts_salmon), group = samplesgroup)

# Filter for genes with low counts across conditions
keep <- rowSums(cpm(fc.dge) > 1) >= 2
fc.dge <- fc.dge[keep, , keep.lib.sizes=FALSE] 

# Normalize for library compositional bias
# calcNormFactors normalizes for RNA composition by finding a set of scaling factors for the library sizes 
# that minimize the log-fold changes between the samples for most genes. 
fc.dge.norm <- calcNormFactors(fc.dge)

# Compute log2 counts-per-million without an offset:
logCPM <- cpm(fc.dge.norm, log = TRUE)
logCPMwonorm <- cpm(fc.dge, log = TRUE)

# Save logCPM to do heatmap:
write.table(logCPM, file = "logCPM_salmon_A_bovis_only.txt", quote = F)


#####################################
# PCA plots
#####################################

metadata2 <- fc.dge.norm$samples
metadata2$batch <- as.character(rep(1,6))
metadata2$gsamples <- rep(c("Bwt", "Bwt", "Bpnca", "Bpnca"), 3)
metadata2$pza <- rep(c("0", "50"), 6)

# By samples - Not normalized
p <- pca(logCPMwonorm, metadata = metadata2)
png(filename = "Pictures/PCAplot_samples_notnorm.png", width = 1600, height = 1200, res = 300)
biplot(p, lab = metadata2$group, labSize = 3, hline = 0, vline = 0, legendPosition = 'right',
       drawConnectors = F, title = "PCA bi-plot", subtitle = "Not normalized, colored by samples",
       colby = 'gsamples', pointSize = 1)
dev.off()

# By batches- not normalized
#png(filename = "Pictures/PCAplot_batches_notnorm.png", width = 1600, height = 1200, res = 300)
#biplot(p, lab = metadata2$group, labSize = 3, hline = 0, vline = 0, legendPosition = 'right',
#       drawConnectors = F, title = "PCA bi-plot", subtitle = "Not normalized, colored by batches",
#       colby = 'batch', pointSize = 1)
#dev.off()


# By samples - Normalized
p <- pca(logCPM, metadata = metadata2)
png(filename = "Pictures/PCAplot_samples_norm.png", width = 1600, height = 1200, res = 300)
biplot(p, lab = metadata2$group, labSize = 3, hline = 0, vline = 0, legendPosition = 'right',
       drawConnectors = F, title = "PCA bi-plot", subtitle = "Normalized and colored by samples",
       colby = 'gsamples', pointSize = 1)
dev.off()

# By batches - Normalized
#png(filename = "Pictures/PCAplot_batches_norm.png", width = 1600, height = 1200, res = 300)
#biplot(p, lab = metadata2$samples, labSize = 3, hline = 0, vline = 0, legendPosition = 'right',
#       drawConnectors = F, title = "PCA bi-plot", subtitle = "Normalized, colored by batches",
#       colby = 'batch', pointSize = 1)
#dev.off()

# By PZA - Normalized
png(filename = "Pictures/PCAplot_pza_norm.png", width = 1600, height = 1200, res = 300)
biplot(p, lab = metadata2$samples, labSize = 3, hline = 0, vline = 0, legendPosition = 'right',
       drawConnectors = F, title = "PCA bi-plot", subtitle = "Normalized, colored by [PZA]",
       colby = 'pza', pointSize = 1)
dev.off()


# Correct batch effect by including batch in model
#Batch <- factor(c("run01", "run01", "run02",
#                  "run01", "run01", "run02",
#                  "run01", "run01", "run02",
#                  "run01", "run01", "run02",
#                  "run02", "run02", "run01",
#                  "run02", "run02", "run01",
#                  "run02", "run02", "run01",
#                  "run02", "run02", "run01",
#                  "run02", "run02", "run03",
#                  "run02", "run02", "run03",
#                  "run02", "run02", "run03",
#                  "run02", "run02", "run03"))

Design <- model.matrix(~0+Treat2)
colnames(Design) <- levels(Treat2)

# Estimate dispersion
fc.dge.disp <- estimateDisp(fc.dge.norm, design = Design, robust = TRUE)
fc.dge.disp$common.dispersion

# MDS plot with color by pza
png(filename = "Pictures/MDSplot_pairwise.png", width = 2732, height = 2048, res = 300)
plotMDS(fc.dge.disp, labels = metadata2$gsamples, 
        col = metadata2$pza)
dev.off()

# MDS plot - basically a PCA
png(filename = "Pictures/MDSplot_common_PCAlike.png", width = 2732, height = 2048, res = 300)
plotMDS(fc.dge.disp, method = "logFC", gene.selection = "common", labels = metadata2$gsamples,
        col = metadata2$pza)
dev.off()

# BCV plot
png(filename = "Pictures/BCVplot.png", width = 1600, height = 1200, res = 300)
plotBCV(fc.dge.disp)
dev.off()

# Fit model / test for DE genes
fit <- glmQLFit(fc.dge.disp, Design, robust = TRUE)
png(filename = "Pictures/QLplot.png", width = 1600, height = 1200, res = 300)
plotQLDisp(fit)
dev.off()


### DE Analysis
# Check if there was a genuine need to adjust for batches by testing for DE between batches
# by testing for DE between all treatments.
#qlf <- glmQLFTest(fit, coef = 2:14)
#tt <- topTags(qlf)
#FDR <- p.adjust(qlf$table$PValue, method="BH")
#sFDR <- sum(FDR < 0.05)
#sdecidetest <- summary(decideTests(qlf))


# Print to txt useful info
########################################
writeLines(c("Common Dispersion:", fc.dge.disp$common.dispersion, "\n", "Toal DE genes at 5% FDR:") , "salmon_DGE_results_A_bovis.txt")


# DE Testing by conditions (comparisons)
########################################
conditions <- read_excel("comparisons.xlsx")

for (i in (1:6)){
  
  file_name  <- paste0(i,"_",conditions[i,1])
  title_name <- as.character(conditions[i,2])
  
  n <- length(conditions[4,])
  x <- as.numeric(unlist(c(conditions[i,3:n])))
  
  
  ###DE Test
  qlf <- glmQLFTest(fit, contrast = x)
  qlfgenes <- rownames(qlf$table)
  qlfx <- qlf$table$logFC
  tt <- topTags(qlf, sort.by = "PValue", n = 10000)
  res <- tt$table
  res <- cbind(res, "FoldChange" = 2^(res$logFC))
  res <- subset(res, select = c(1,2,7,3,4,5,6))
  write.xlsx(res, paste0("Tables/",file_name,".xlsx"), 
             col.names=TRUE, row.names=FALSE, append=FALSE)
  if (nrow(res[res$FDR < 0.05,]) > 0) {
    write.xlsx(res[res$FDR<0.05,],paste0("Tables/FDR/",file_name,".xlsx"))
  }
  cat("\n", file = "salmon_DGE_results_A_bovis.txt", append = TRUE)
  write.table(summary(decideTests(qlf)), file = "salmon_DGE_results_A_bovis.txt", append = T)
  
  ###Plots
  ##MD
  png(filename = paste0("Pictures/",file_name,".png"))
  plotMD(qlf, main = title_name)
  abline(h = c(-1,1), col = "blue")
  dev.off()
  
  
  ##Volcano Plot A
  png(filename = paste0("Pictures/VolcanoA/",file_name,"_VolcanoA.png"), width = 2732, height = 2048, res = 300)
  with(res, plot(res$logFC, -log10(FDR), pch=20, main=title_name,
                 xlim=c(-13, 6), ylim = c(0.5, 3.5), xlab = "log2FC"))
  abline(v = c(-1.5, 1.5), col = "blue")
  abline(h = c(1, 1.5), col = "blue")
  with(subset(res, logFC >= 1.5 & -log10(FDR) >= 1), points(logFC, -log10(FDR), pch=20, col="darkred"))
  with(subset(res, logFC <= -1.5 & -log10(FDR) >= 1), points(logFC, -log10(FDR), pch=20, col="forestgreen"))
  with(subset(res, logFC >= 1.5 & -log10(FDR) >= 1.5), points(logFC, -log10(FDR), pch=20, col="red"))
  with(subset(res, logFC <= -1.5 & -log10(FDR) >= 1.5), points(logFC, -log10(FDR), pch=20, col="green"))
  with(subset(res, logFC >= 1.5 & -log10(FDR) >= 1.5), textxy(logFC, -log10(FDR), labs=genes, cex=.5))
  with(subset(res, logFC <= -1.5 & -log10(FDR) >= 1.5), textxy(logFC, -log10(FDR), labs=genes, cex=.5))
  dev.off()
  GIDs <- c(subset(res$genes, res$logFC > 1.5 & -log10(res$FDR) > 1),
            subset(res$genes, res$logFC < -1.5 & -log10(res$FDR) > 1))
  cat(c("\n",title_name," gene names (Volcano Plot A):", "\n", GIDs),
      file =  "salmon_DGE_results_A_bovis.txt", append = T)
  
  ##Volcano Plot B
  png(filename = paste0("Pictures/VolcanoB/",file_name,"_VolcanoB.png"), width = 2732, height = 2048, res = 300)
  with(res, plot(res$logFC, -log10(FDR), pch=20, main=title_name,
                 xlim=c(-13, 6), ylim = c(0.5, 3.5), xlab = "log2FC"))
  abline(v = c(-2, 2), col = "blue")
  abline(h = c(1, 1.5), col = "blue")
  with(subset(res, logFC >= 2 & -log10(FDR) >= 1), points(logFC, -log10(FDR), pch=20, col="darkred"))
  with(subset(res, logFC <= -2 & -log10(FDR) >= 1), points(logFC, -log10(FDR), pch=20, col="forestgreen"))
  with(subset(res, logFC >= 2 & -log10(FDR) >= 1.5), points(logFC, -log10(FDR), pch=20, col="red"))
  with(subset(res, logFC <= -2 & -log10(FDR) >= 1.5), points(logFC, -log10(FDR), pch=20, col="green"))
  with(subset(res, logFC >= 2 & -log10(FDR) >= 1.5), textxy(logFC, -log10(FDR), labs=genes, cex=.5))
  with(subset(res, logFC <= -2 & -log10(FDR) >= 1.5), textxy(logFC, -log10(FDR), labs=genes, cex=.5))
  dev.off()
  GIDs <- c(subset(res$genes, res$logFC > 2 & -log10(res$FDR) > 1),
            subset(res$genes, res$logFC < -2 & -log10(res$FDR) > 1))
  cat(c("\n",title_name," gene names (Volcano Plot A):", "\n", GIDs, "\n\n"), 
      file =  "salmon_DGE_results_A_bovis.txt", append = T)
  
  
  # Save data for heatmaps
  assign(paste("tgenesFDRA", i, sep = "_"), subset(res$genes, res$FDR < 0.1))
  assign(paste("tgenesFDRB", i, sep = "_"), subset(res$genes, res$FDR < 0.05))
  assign(paste("tgenesPV", i, sep = "_"), subset(res$genes, res$PValue < 0.01))
  assign(paste("qlf", i, sep = ""), qlfx)
  
  #  write.table(subset(res$genes, res$FDR < 0.1), file = paste0("heatmapfiles/",i,"_FDRA_",file_name,".txt"), sep=",", col.names = FALSE)
  #  write.table(subset(res$genes, res$FDR < 0.05), file = paste0("heatmapfiles/",i,"_FDRB_",file_name,".txt"), sep=",", col.names = FALSE)
  #  write.table(subset(res$genes, res$PValue < 0.01), file = paste0("heatmapfiles/",i,"_PV2_",file_name,".txt"), sep=",", col.names = FALSE)                        
}            


# Save data for heatmaps
########################################

# Make a tab-delimited table with logFC values per run
qlfgenes <- rownames(qlf$table)
df2 <- data.frame("Genes" = qlfgenes, "q1" = qlf1, "q2" = qlf2, "q3" = qlf3, "q4" = qlf4, "q5" = qlf5, "q6" = qlf6)
write.table(df2, file = "salmon_logFC.txt")


# Make a tab-delimited table with names of top genes per run
# Pval < 0.01
TG_PV_all <- list("q1" = tgenesPV_1, "q2" = tgenesPV_2, "q3" = tgenesPV_3, "q4" = tgenesPV_4, "q5" = tgenesPV_5, "q6" = tgenesPV_6)
# FDR A: FDR < 0.1
TG_FDRa_all <- list("q1" = tgenesFDRA_1, "q2" = tgenesFDRA_2, "q3" = tgenesFDRA_3, "q4" = tgenesFDRA_4, "q5" = tgenesFDRA_5, "q6" = tgenesFDRA_6)
# FDR B: FDR < 0.05
TG_FDRb_all <- list("q1" = tgenesFDRB_1, "q2" = tgenesFDRB_2, "q3" = tgenesFDRB_3, "q4" = tgenesFDRB_4, "q5" = tgenesFDRB_5, "q6" = tgenesFDRB_6)


save.image("Adiana_bovis_only.RData")










