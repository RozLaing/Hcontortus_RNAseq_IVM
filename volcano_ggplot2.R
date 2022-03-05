# R code to make volcano plots with ggplot2

# input file is summary of all DESEQ2 results in format:
# M_MHco18vMHco3     baseMean log2FoldChange      lfcSE       stat      pvalue        padj ...
# HCON_00000010    0.7135306     -2.6366072 3.90184133 -0.6757341 0.499209485          NA ...
# HCON_00000020 2841.1333170      0.2289005 0.06252609  3.6608796 0.000251351 0.001240881 ...
# HCON_00000030 1271.9827950      0.3185192 0.10575708  3.0118002 0.002597034 0.009507994 ...
# ...


library(ggplot2)
DE<-read.table("summary_allmales_allfemales.txt", header=T, sep="\t", na.strings="NA")

# males
SIG <- subset(DE, padj<0.01 & padj.1>0.01 & padj.2>0.01)
IVM1 <- subset(DE, padj<0.01 & padj.1>0.01 & padj.2<0.01)
IVM2 <- subset(DE, padj<0.01 & padj.1<0.01 & padj.2<0.01)

ggplot(SIG, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(alpha = 0.5, col = "#0075DC") +  scale_colour_manual(values=c( "#F0A3FF", "gold")) +
  geom_point(data=IVM1, aes(x = log2FoldChange, y = (-log10(pvalue)), color='IVM1', alpha = 1)) +
  geom_point(data=IVM2, aes(x = log2FoldChange, y = (-log10(pvalue)), color='IVM2', alpha = 1)) +
  guides(alpha = FALSE, color = F) +
  labs(x = "log 2 fold change MHco18vMHco3 males", y = "-log10 pvalue") + ylim(0,180) + xlim(-12,12) 


# females
SIGF <- subset(DE, padj.5<0.01 & padj.7>0.01)
IVMF <- subset(DE, padj.5<0.01 & padj.7<0.01)

ggplot(SIGF, aes(x = log2FoldChange.5, y = -log10(pvalue.5))) + 
  geom_point(alpha = 0.5, col = "#0075DC") +  scale_colour_manual(values=c( "#F0A3FF")) +
  geom_point(data=IVMF, aes(x = log2FoldChange.5, y = (-log10(pvalue.5)), color='IVM1', alpha = 1)) +
  guides(alpha = FALSE, color = F) +
  labs(x = "log 2 fold change MHco18vMHco3 females", y = "-log10 pvalue") +  ylim(0,180) + xlim(-12,12) 


