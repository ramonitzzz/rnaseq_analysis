#setwd("~/OneDrive/Bioinfo/Proyectos/Tuberculosis/RNA-seq/Run5/Fernando")
setwd("~/Desktop/20_RNAseq/Adiana/H37Rv_only")



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

#dir.create("heatmapfiles")
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
samplesinfo <- c("IH-p0;pH-6-3_PZA-0", "IH-p50;pH-6-3_PZA-50",
                 "IIH-p0;pH-6-3_PZA-0", "IIH-p50;pH-6-3_PZA-50",
                 "IIIH-p0;pH-6-3_PZA-0", "IIIH-p50;pH-6-3_PZA-50")

# Try to group by treat + sample

Treat2 <- factor(c("H37Rv_pH-6-3_PZA-0", "H37Rv_pH-6-3_PZA-50",
                   "H37Rv_pH-6-3_PZA-0", "H37Rv_pH-6-3_PZA-50",
                   "H37Rv_pH-6-3_PZA-0", "H37Rv_pH-6-3_PZA-50"))


samplesgroup <-c("H1", "H1", "H2", "H2", "H3", "H3")


# Counts
counts_salmon = read.table("counts_h37rv.tsv", sep="\t", header = TRUE, row.names = 1, stringsAsFactors = F)

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
write.table(logCPM, file = "logCPM_salmon_A_H37Rv_only.txt", quote = F)


#####################################
# PCA plots
#####################################

metadata2 <- fc.dge.norm$samples
metadata2$batch <- as.character(rep(1,6))
metadata2$gsamples <- rep(c("H"), each=6)
metadata2$pza <- rep(c("0", "50"), 3)

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

### Prepare objects for GOseq analysis

table_length <- read.table("genes_length_H37Rv.tsv", sep="\t", col.names=c("Gene","Length"), stringsAsFactors = F)
Size_genes <- as.integer(table_length$Length)
names(Size_genes) <- table_length$Gene

# Variable GOterm containing GO terms 
GOT <- read.table("allGoTermsForMycobacterium_tuberculosis_H37Rv.txt", sep = "\t", stringsAsFactors = F)
GOTnames <- read.table("names.tsv", sep = "\t", stringsAsFactors = F)
GOTnames <- GOTnames$V1
GOT <- GOT[-1,]
i <- 1
X <- c()
Y <- c()
Name <- c()
GOterm <- list()
while (i < length(GOT$V1)) {
  X <- c(X, GOT$V1[i])
  Y <- c(Y, paste0("GO:", GOT$V2[i]))
  
  if (GOT$V1[i] != GOT$V1[i+1] & GOT$V4[i] == "biological_process") {
    #  if (GOT$V1[i] != GOT$V1[i+1]) {
    W <- list("X" = Y)
    Name <- c(Name, GOT$V1[i+1])
    X <- c()
    Y <- c()
    GOterm  <- c(GOterm, W)
  }
  i <- i + 1
}
names(GOterm) <- Name




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
writeLines(c("Common Dispersion:", fc.dge.disp$common.dispersion, "\n", "Toal DE genes at 5% FDR:") , "salmon_DGE_results_A_H37Rv.txt")


# DE Testing by conditions (comparisons)
########################################
file_name <- "H37Rv_pza0vs50_pH6-3_FDR"
title_name <- "H37Rv pH 6.3, PZA 0 vs. 50"

qlf <- glmQLFTest(fit, contrast = c(-1,1))
qlfgenes <- rownames(qlf$table)
qlf1 <- qlf$table$logFC
tt <- topTags(qlf, sort.by = "PValue", n = 10000)
res <- tt$table
res <- cbind(res, "FoldChange" = 2^(res$logFC))
write.xlsx(res, "Tables/H37Rv_pza0vs50_pH6-3.xlsx", col.names=TRUE, row.names=FALSE, append=FALSE)
if (nrow(res[res$FDR < 0.05,]) > 0) {
  write.xlsx(res[res$FDR<0.05,], "H37Rv_pza0vs50_pH6-3_FDR.xlsx")
}

write.table(summary(decideTests(qlf)), file = "salmon_DGE_results_A_H37Rv.txt", append = T)

###Gene enrichment analysis (GEA) - GOseq (based on pValue)
#Table with genes (column 1: genes; column 2: 1(significant) or 0 (non-significant))
selected <- c(rep(1,nrow(res[res$PValue<=0.05,])), rep(0,nrow(res[res$PValue>0.05,])))
DE <- matrix(selected, ncol=1)
row.names(DE) <- res$genes

Size <- Size_genes[names(Size_genes)%in%rownames(DE)]
DE <- DE[rownames(DE)%in%names(Size_genes)]

Size_genes <- Size_genes[names(Size_genes)%in%names(GOterm)]
GOterm <- GOterm[names(GOterm)%in%names(Size_genes)]

Size <- Size_genes[names(Size_genes)%in%res$genes]
DE <- DE[res$genes%in%names(Size_genes)]
names(DE) <- res$genes[res$genes%in%names(Size_genes)]

DE[is.na(DE)] = 0

pwf = nullp(DE, bias.data=Size)
GO.wall=goseq(pwf,gene2cat=GOterm,use_genes_without_cat=FALSE)
dev.off()

#Analysis by Biological proc. (BP), Molecular funct. (MF) or Cellular compartment (CC)

#Over-represented pathways (Biological proc. bp)
GO.wall.or.bp    <- GO.wall[GO.wall$ontology == "BP" & GO.wall$over_represented_pvalue < 0.05,]
GO.wall.or.bp    <- na.omit(GO.wall.or.bp)
GO.wall.or.bp.df <- data.frame("Description"            = as.character(GO.wall.or.bp$term),
                               "Number of DE Genes"     = as.integer(GO.wall.or.bp$numDEInCat),
                               "Number of Genes in Cat" = as.integer(GO.wall.or.bp$numInCat))


GO.wall.or.bp.df <- data.frame(cbind(GO.wall.or.bp.df, 
                                     "Gene ratio" = round(GO.wall.or.bp.df$Number.of.DE.Genes/GO.wall.or.bp.df$Number.of.Genes.in.Cat, digits = 2)),
                               "PValue" = GO.wall.or.bp$over_represented_pvalue)
GO.wall.or.bp.df <- na.omit(GO.wall.or.bp.df)

if (nrow(GO.wall.or.bp.df) > 0) {
  write.xlsx(GO.wall.or.bp.df, paste0("Tables/GOseq/bp/",file_name,".xlsx"), row.names=FALSE)
  
  
  
  
  ##Plot GEA
  #Dot plot
  tiff(filename = paste0("Pictures/GEA/bp/",file_name,".tiff"), width=3024, height= 4048, res=300)
  print (ggplot(GO.wall.or.bp.df) +
           geom_point(mapping = aes(x=Gene.ratio, y = reorder(Description, Gene.ratio), color = PValue)) + 
           scale_color_gradientn(colors=c("#FC4E07", "#E7B800", "#00AFBB")) +
           labs(x = "Gene ratio", y = "" ) +
           ggtitle("Over-representated GO terms") +
           theme_minimal() +
           theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  #Pie chart
  pal<-c(randomColor(count=nrow(GO.wall.or.bp.df)))
  tiff(filename = paste0("Pictures/GEA/bp/piechart/",file_name,".tiff"), width=6024, height= 2048, res=300)
  print (ggplot(GO.wall.or.bp.df, aes (x="", y=Number.of.DE.Genes, fill = factor(Description)))+
           geom_col(position = 'stack', width = 0.1)  +
           #    scale_color_gradientn(colors=c("#FC4E07", "#E7B800", "#00AFBB")) + 
           scale_fill_manual(values=pal) +
           coord_polar("y") +
           theme_classic() +
           theme(plot.title = element_text(hjust = 0.5),
                 axis.line  = element_blank(),
                 axis.text  = element_blank(),
                 axis.ticks = element_blank(),
                 legend.position = "right") +
           labs(fill="", x = NULL, y = NULL))
  dev.off()
}


#Over-represented pathways (Molecular func. mf)
GO.wall.or.mf    <- GO.wall[GO.wall$ontology == "MF" & GO.wall$over_represented_pvalue < 0.05,]
GO.wall.or.mf    <- na.omit(GO.wall.or.mf)
GO.wall.or.mf.df <- data.frame("Description"            = as.character(GO.wall.or.mf$term),
                               "Number of DE Genes"     = as.integer(GO.wall.or.mf$numDEInCat),
                               "Number of Genes in Cat" = as.integer(GO.wall.or.mf$numInCat))


GO.wall.or.mf.df <- data.frame(cbind(GO.wall.or.mf.df, 
                                     "Gene ratio" = round(GO.wall.or.mf.df$Number.of.DE.Genes/GO.wall.or.mf.df$Number.of.Genes.in.Cat, digits = 2)),
                               "PValue" = GO.wall.or.mf$over_represented_pvalue)
GO.wall.or.mf.df <- na.omit(GO.wall.or.mf.df)

if (nrow(GO.wall.or.mf.df)) {
  write.xlsx(GO.wall.or.mf.df, paste0("Tables/GOseq/mf/",file_name,"_bp.xlsx"), row.names=FALSE) 
  
  
  
  ##Plot GEA
  #Dot plot
  tiff(filename = paste0("Pictures/GEA/mf/",file_name,".tiff"), width=3024, height= 4048, res=300)
  print (ggplot(GO.wall.or.mf.df) +
           geom_point(mapping = aes(x=Gene.ratio, y = reorder(Description, Gene.ratio), color = PValue)) + 
           scale_color_gradientn(colors=c("#FC4E07", "#E7B800", "#00AFBB")) +
           labs(x = "Gene ratio", y = "" ) +
           ggtitle("Over-representated GO terms") +
           theme_minimal() +
           theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  #Pie chart
  pal<-c(randomColor(count=nrow(GO.wall.or.mf.df)))
  tiff(filename = paste0("Pictures/GEA/mf/piechart/",file_name,".tiff"), width=6024, height= 2048, res=300)
  print (ggplot(GO.wall.or.mf.df, aes (x="", y=Number.of.DE.Genes, fill = factor(Description)))+
           geom_col(position = 'stack', width = 0.1)  +
           #    scale_color_gradientn(colors=c("#FC4E07", "#E7B800", "#00AFBB")) + 
           scale_fill_manual(values=pal) +
           coord_polar("y") +
           theme_classic() +
           theme(plot.title = element_text(hjust = 0.5),
                 axis.line  = element_blank(),
                 axis.text  = element_blank(),
                 axis.ticks = element_blank(),
                 legend.position = "right") +
           labs(fill="", x = NULL, y = NULL))
  dev.off()
}

#Over-represented pathways (Molecular func. mf)
GO.wall.or.cc    <- GO.wall[GO.wall$ontology == "CC" & GO.wall$over_represented_pvalue < 0.05,]
GO.wall.or.cc    <- na.omit(GO.wall.or.cc)
GO.wall.or.cc.df <- data.frame("Description"            = as.character(GO.wall.or.cc$term),
                               "Number of DE Genes"     = as.integer(GO.wall.or.cc$numDEInCat),
                               "Number of Genes in Cat" = as.integer(GO.wall.or.cc$numInCat))


GO.wall.or.cc.df <- data.frame(cbind(GO.wall.or.cc.df, 
                                     "Gene ratio" = round(GO.wall.or.cc.df$Number.of.DE.Genes/GO.wall.or.cc.df$Number.of.Genes.in.Cat, digits = 2)),
                               "PValue" = GO.wall.or.cc$over_represented_pvalue)
GO.wall.or.cc.df <- na.omit(GO.wall.or.cc.df)


if (nrow(GO.wall.or.cc.df)){
  tryCatch(write.xlsx(GO.wall.or.cc.df, paste0("Tables/GOseq/cc/",file_name,"_bp.xlsx"), row.names=FALSE)) 
  ##Plot GEA
  #Dot plot
  tiff(filename = paste0("Pictures/GEA/cc/",file_name,".tiff"), width=3024, height= 4048, res=300)
  print (ggplot(GO.wall.or.cc.df) +
           geom_point(mapping = aes(x=Gene.ratio, y = reorder(Description, Gene.ratio), color = PValue)) + 
           scale_color_gradientn(colors=c("#FC4E07", "#E7B800", "#00AFBB")) +
           labs(x = "Gene ratio", y = "" ) +
           ggtitle("Over-representated GO terms") +
           theme_minimal() +
           theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  #Pie chart
  pal<-c(randomColor(count=nrow(GO.wall.or.cc.df)))
  tiff(filename = paste0("Pictures/GEA/cc/piechart/",file_name,".tiff"), width=6024, height= 2048, res=300)
  print (ggplot(GO.wall.or.cc.df, aes (x="", y=Number.of.DE.Genes, fill = factor(Description)))+
           geom_col(position = 'stack', width = 0.1)  +
           #    scale_color_gradientn(colors=c("#FC4E07", "#E7B800", "#00AFBB")) + 
           scale_fill_manual(values=pal) +
           coord_polar("y") +
           theme_classic() +
           theme(plot.title = element_text(hjust = 0.5),
                 axis.line  = element_blank(),
                 axis.text  = element_blank(),
                 axis.ticks = element_blank(),
                 legend.position = "right") +
           labs(fill="", x = NULL, y = NULL))
  dev.off()
}



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
    file =  "salmon_DGE_results_A_H37Rv.txt", append = T)

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
    file =  "salmon_DGE_results_A_H37Rv.txt", append = T)

write.table(subset(res$genes, res$FDR < 0.1), file = paste0("heatmapfiles/FDRA_",file_name,".txt"), sep=",", col.names = FALSE)
write.table(subset(res$genes, res$FDR < 0.05), file = paste0("heatmapfiles/FDRB_",file_name,".txt"), sep=",", col.names = FALSE)
write.table(subset(res$genes, res$PValue < 0.01), file = paste0("heatmapfiles/PV2_",file_name,".txt"), sep=",", col.names = FALSE)                        
           


save.image("Adiana_H37Rv_only.RData")










