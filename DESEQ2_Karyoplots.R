# R code to perform DESEQ2 analysis, make MA and PCA plots, heatmaps, plot normalised counts and karyoplots

# input files are:
# featureCounts.tabular in format:
# MHco18_1 MHco18_2 MHco18_3 MHco3_1 MHco3_2 MHco3_3 F2IVM_1 F2IVM_2 F2IVM_3 ...
# HCON_00000010        0        0        0       2       0       2       0       0       4 ...
# HCON_00000020     2692     2300     2538    2524    2870    2396    2925    3099    2245 ...
# HCON_00000030     1168      966      969     920    1127     972    1508    1622    1090 ...
# ...

# metadata.txt in format:
# Strain Sex Trizol
# MHco18_1 MHco18   M      A
# MHco18_2 MHco18   M      B
# MHco18_3 MHco18   M      C
# MHco3_1   MHco3   M      A
# MHco3_2   MHco3   M      B
# MHco3_3   MHco3   M      C
#...

# DESEQ2 with batch correction (provided by Bruce Rosa), alpha set to 0.01
library("DESeq2")
library("geneplotter")
COUNTS <- as.matrix(read.table(file="featureCounts_allmales_WBPS15gff.tabular", sep="\t", header=TRUE, row.names=1))
META <- read.table(file="allmales_metadata.txt", sep="\t", header=TRUE, row.names=1)
dds <- DESeqDataSetFromMatrix(countData=COUNTS, colData=META, design = ~Strain)
# for batch correction of female samples: dds <- DESeqDataSetFromMatrix(countData=COUNTS, colData=META, design = ~Trizol + Strain)
dds <- DESeq(dds)
outputtable <- results(dds, alpha = 0.01, contrast=c("Strain", "F2IVM", "F2CTL"))
summary (outputtable)
write.table(outputtable, file="DESEQ_WBPS15gff_F2IVMvF2CTLmales_alpha0.01.tabular", sep="\t", col.names = NA)

# MA plot - this is based on the specified pairwise comparison from DESEQ2 results
# LEGEND = The log2 fold change for a particular comparison is plotted on the y-axis and the average of the counts normalized by size factor is shown on the x-axis. Each gene is represented with a dot. Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.
res <- results(dds, alpha = 0.01, contrast=c("Strain", "F2IVM", "F2CTL"))
plotMA(res)
plotMA(res, ylim=c(-12,12))
#Dispersion estimates
plotDispEsts(dds)
mcols(res, use.names = TRUE)

# PCA - DESeq-based PCA - used top 500 variable genes
# LEGEND = PCA plot with rlog transformed data
library(ggplot2)
rld <- rlogTransformation(dds, blind=TRUE)
data <- plotPCA(rld, intgroup=c("Strain", "Trizol"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=Strain, shape=Trizol)) +
  geom_point(size=3) + scale_colour_manual(values=c( "#F0A3FF", "gold", "#0075DC", "turquoise4", "orange")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + xlim(-20,20) + ylim(-20,20)

# heatmap - adapted from Bioconductor RNAseq workflow (Love, Anders. Kim and Huber, October 2019: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)
# LEGEND = Heatmap showing the Euclidean distances between the samples as calculated from the regularized log transformation.
library(pheatmap)
library("RColorBrewer")
# transpose rld and obtain Euclidean distance
# use assay accessor to obtain rlogged observations
dist_rl = dist(t(assay(rld)))
#convert to a matrix
dist_mat = as.matrix(dist_rl)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(dist_mat,
         clustering_distance_rows = dist_rl,
         clustering_distance_cols = dist_rl,
         col = colors)

# for comparison, heatmap using Poisson Distance (measure of dissimilarity between counts, also takes into account inherent variance of counts)
library("PoiClaClu")
poisd <- PoissonDistance(t(COUNTS))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- c("MHco18_1", "MHco18_2", "MHco18_3", "MHco3_1", "MHco3_2", "MHco3_3", "F2IVM_1", "F2IVM_2", "F2IVM_3", "F2CTL1", "F2CTL_2", "F2CTL_3", "F2BZ_1", "F2BZ_2", "F2BZ_3")
colnames(samplePoisDistMatrix) <- c("MHco18_1", "MHco18_2", "MHco18_3", "MHco3_1", "MHco3_2", "MHco3_3", "F2IVM_1", "F2IVM_2", "F2IVM_3", "F2CTL1", "F2CTL_2", "F2CTL_3", "F2BZ_1", "F2BZ_2", "F2BZ_3")
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# counts plot from RNAseq workflow (Love, Anders. Kim and Huber, October 2019: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)
# LEGEND = Normalized counts for a single gene over treatment group
# plain version
plotCounts(dds, gene = "HCON_00162780", intgroup=c("Strain"))
# pretty version
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = "HCON_00162780", intgroup=c("Strain"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Strain, y = count, color = Strain)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

# Karyoplots - based on KaryploteR tutorial (http://bernatgel.github.io.convey.pro/l/2PzJRja)
library(karyoploteR)

# make custom genome for Hc based on chromosome lengths (for ChV-only plots, just remove other chromosomes in custom genome and proceed as for full analysis)
custom.genome <- toGRanges(data.frame(chr=c("chr1", "chr2", "chr3", "chr4", "chr5", "chrX"), start=c(1,1,1,1,1,1), end=c(45774638,47382675,43560351,51826578,48868367,46012676)))
#check looks ok
kp <- plotKaryotype(plot.type=1, custom.genome,labels.plotter=NULL) 
kpAddChromosomeNames(kp, chr.names=c("chI", "chII", "chIII", "chIV", "chV", "chX"))

# convert gff3 file to txdb file
library(GenomicFeatures)
Hc_txdb <- makeTxDbFromGFF("~/Documents/2020/RNAseq_NewGFF_April2020/HCON_V4_curated_20200422_WBPS15_cleannames.gff3", format="gff3")
Hc_txdb

# extract data for genes
Hc.genes <- genes(Hc_txdb)
Hc.genes

# make sure gene names in DESEQ outputtable match gff3 (e.g. remove gene:)
res <- outputtable
res

mcols(Hc.genes) <- res[names(Hc.genes), c("log2FoldChange", "stat", "pvalue", "padj")]
head(Hc.genes, n=4)

# to plot P values of all sig genes
# first filter out genes with NA in Padj column
filtered.Hc.genes <- Hc.genes[!is.na(Hc.genes$padj)]
log.pval <- -log10(filtered.Hc.genes$padj)
mcols(filtered.Hc.genes)$log.pval <- log.pval
filtered.Hc.genes

# filter for sig genes
sign.genes <- filtered.Hc.genes[filtered.Hc.genes$padj < 0.01,]

# filter on log2 fold change
range(sign.genes$log2FoldChange)
fc.ymax <-ceiling(max(abs(range(sign.genes$log2FoldChange))))
fc.ymin <- -fc.ymax

# size points to represent p value
cex.val <- sqrt(sign.genes$log.pval)/3

# colour points
col.over <- "#FFBD07AA"
col.under <- "#00A6EDAA"
sign.col <- rep(col.over, length(sign.genes))
sign.col[sign.genes$log2FoldChange<0] <- col.under

# basic plot 
kp <- plotKaryotype(custom.genome,labels.plotter=NULL, plot.type=1)
kpAddChromosomeNames(kp, chr.names=c("chI", "chII", "chIII", "chIV", "chV", "chX"))
kpAddBaseNumbers(kp)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, col=sign.col)
#add axis and labels
kpAxis(kp, ymax=fc.ymax, ymin = fc.ymin)
kpAddLabels(kp,labels = "log2 FC", srt=90, pos=1, label.margin = 0.04)

# multiple plots with FST panel

# additional input files are:
# FST data (from Stephen Doyle) in format:
# CHR	START	END	FST_R1
# chr1	2500	7500	0.00415596
# chr1	7500	12500	0.00886197
# ...

# data=sign.genes for each pairwise comparison as per basic plot above

# convert FST data
Fdata <- read.table("FST_karyoR1.txt", header = T)
Fdata <- makeGRangesFromDataFrame(Fdata, keep.extra.columns=TRUE)
Fdata <- toGRanges(Fdata, genome = custom.genome)

# plot data as panels (r0 and r1)
kp <- plotKaryotype(custom.genome,labels.plotter=NULL, plot.type = 4) 
kpAddChromosomeNames(kp, chr.names=c("chI", "chII", "chIII", "chIV", "chV", "chX"))

kpDataBackground(kp, r0=0, r1=0.12)
kp <- kpPoints(kp, data=Fdata, y=Fdata$FST_R1, ymin=0, ymax=0.08, col="gray62", r0=0, r1=0.12)
kpAxis(kp, ymax=0.08, ymin=0, cex = 0.8, numticks = 2, r0=0, r1=0.12)

kpDataBackground(kp, r0=0.15, r1=0.29)
kpPoints(kp, data=M_F2IVMvCTL, y=M_F2IVMvCTL$log2FoldChange, cex=M_F2IVMvCTL.cex.val, ymax=12, ymin=-12, col=M_F2IVMvCTL.sign.col, r0=0.15, r1=0.29)
kpAxis(kp, ymax=12, ymin=-12, cex = 0.8, r0=0.15, r1=0.29)

kpDataBackground(kp, r0=0.32, r1=0.46)
kpPoints(kp, data=F_F2IVMv3, y=F_F2IVMv3$log2FoldChange, cex=F_F2IVMv3.cex.val, ymax=12, ymin=-12, col=F_F2IVMv3.sign.col, r0=0.32, r1=0.46)
kpAxis(kp, ymax=12, ymin=-12, cex = 0.8, r0=0.32, r1=0.46)

kpDataBackground(kp, r0=0.49, r1=0.63)
kpPoints(kp, data=M_F2IVMv3, y=M_F2IVMv3$log2FoldChange, cex=M_F2IVMv3.cex.val, ymax=12, ymin=-12, col=M_F2IVMv3.sign.col, r0=0.49, r1=0.63)
kpAxis(kp, ymax=12, ymin=-12, cex = 0.8, r0=0.49, r1=0.63)

kpDataBackground(kp, r0=0.66, r1=0.80)
kpPoints(kp, data=F_MHco18v3, y=F_MHco18v3$log2FoldChange, cex=F_MHco18v3.cex.val, ymax=12, ymin=-12, col=F_MHco18v3sign.col, r0=0.66, r1=0.80)
kpAxis(kp, ymax=12, ymin=-12, cex = 0.8, r0=0.66, r1=0.80)

kpDataBackground(kp, r0=0.83, r1=0.97)
kpPoints(kp, data=M_Hco18v3, y=M_Hco18v3$log2FoldChange, cex=M_Hco18v3.cex.val, ymax=12, ymin=-12, col=M_Hco18v3.sign.col, r0=0.83, r1=0.97)
kpAxis(kp, ymax=12, ymin=-12, cex = 0.8, r0=0.83, r1=0.97)


