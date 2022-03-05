#R script to generate scatter plots

# input file (summary_allmales_allfemales.txt) is summary of all DESEQ2 results in format:
# M_MHco18vMHco3     baseMean log2FoldChange      lfcSE       stat      pvalue        padj ...
# HCON_00000010    0.7135306     -2.6366072 3.90184133 -0.6757341 0.499209485          NA ...
# HCON_00000020 2841.1333170      0.2289005 0.06252609  3.6608796 0.000251351 0.001240881 ...
# HCON_00000030 1271.9827950      0.3185192 0.10575708  3.0118002 0.002597034 0.009507994 ...
# ...

#input file (allIVMR.txt) is above filtered for sig DE (padj<0.01) and same direction of logfoldchange in all IVMR pairwise comparisons

library(ggplot2)
library(ggrepel)

DE<-read.table("summary_allmales_allfemales.txt", header=T, sep="\t", na.strings="NA")
all <- as.vector(read.table(file="allIVMR.txt", header = TRUE, na.strings="NA")) 
s1 <- subset(Z, X == "HCON_00155390") #cky-1

#scatterplots for pairwise comparisons
#parental pairwise comparisons - males x females
Z <- subset(DE, padj<0.01 & padj.5>0.01)

P1 <- ggplot(Z, aes(x = log2FoldChange, y = log2FoldChange.5)) + 
  geom_point(alpha = 0.5, col = "gray80") +  scale_colour_manual(values=c("#0075DC")) +
  geom_point(data=all, aes(x = log2FoldChange, y = log2FoldChange.5, color='IVM', alpha = 0.4)) +
  labs(
    x = "log 2 fold change MHco18vMHco3 males",
    y = "log 2 fold change MHco18vMHco3 females") + xlim(-10, 10) + ylim(-10,10) +
  guides(alpha = FALSE) +
  geom_text_repel(data=s1, label = "HCON_00155390", hjust=1,vjust=1, nudge_x = 2, nudge_y = 0.5, cex = 3, col = "black") + geom_point(data=s1, col = "#0075DC") + theme(legend.position = "none")

#parental males x genetic cross males
Y <- subset(DE, padj<0.01 & padj.1>0.01)

P2 <- ggplot(Y, aes(x = log2FoldChange, y = log2FoldChange.1)) + 
  geom_point(alpha = 0.5, col = "gray80") +  scale_colour_manual(values=c("#0075DC")) +
  geom_point(data=all, aes(x = log2FoldChange, y = log2FoldChange.1, color='IVM', alpha = 0.4)) +
  labs(
    x = "log 2 fold change MHco18vMHco3 males",
    y = "log 2 fold change F2IVMvF2CTL males") + xlim(-10, 10) + ylim(-10,10) +
  guides(alpha = FALSE) +
  geom_text_repel(data=s1, label = "HCON_00155390", hjust=1,vjust=1, nudge_x = 2, nudge_y = 0.5, cex = 3, col = "black") + geom_point(data=s1, col = "#0075DC") + theme(legend.position = "none")

#parental males x genetic cross v ISE males
X <- subset(DE, padj<0.01 & padj.2>0.01)

P3 <- ggplot(X, aes(x = log2FoldChange, y = log2FoldChange.2)) + 
  geom_point(alpha = 0.5, col = "gray80") +  scale_colour_manual(values=c("#0075DC")) +
  geom_point(data=all, aes(x = log2FoldChange, y = log2FoldChange.2, color='IVM', alpha = 0.4)) +
  labs(
    x = "log 2 fold change MHco18vMHco3 males",
    y = "log 2 fold change F2IVMvMHco3 males") + xlim(-10, 10) + ylim(-10,10) +
  guides(alpha = FALSE) +
  geom_text_repel(data=s1, label = "HCON_00155390", hjust=1,vjust=1, nudge_x = 2, nudge_y = 0.5, cex = 3, col = "black") + geom_point(data=s1, col = "#0075DC") + theme(legend.position = "none")

#parental females x genetic cross v ISE females
W <- subset(DE, padj.5<0.01 & padj.7>0.01)

P4 <- ggplot(W, aes(x = log2FoldChange.5, y = log2FoldChange.7)) + 
  geom_point(alpha = 0.5, col = "gray80") +  scale_colour_manual(values=c("#0075DC")) +
  geom_point(data=all, aes(x = log2FoldChange.5, y = log2FoldChange.7, color='IVM', alpha = 0.4)) +
  labs(
    x = "log 2 fold change MHco18vMHco3 females",
    y = "log 2 fold change F2IVMvMHco3 females") + xlim(-10, 10) + ylim(-10,10) +
  guides(alpha = FALSE) +
  geom_text_repel(data=s1, label = "HCON_00155390", hjust=1,vjust=1, nudge_x = 2, nudge_y = 0.5, cex = 3, col = "black") + geom_point(data=s1, col = "#0075DC") + theme(legend.position = "none")

library("gridExtra")
grid.arrange(P1, P2, P3, P4, nrow = 2, ncol = 2)

