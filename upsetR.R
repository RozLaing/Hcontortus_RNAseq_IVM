#R code to generate UpSet plots

#input files are lists of differentally expressed genes as .txt files

#script from http://genomespot.blogspot.com/2017/09/upset-plots-as-replacement-to-venn.html
library(plyr)
library(reshape2)
library(UpSetR)

# make a list of input files to be read - save as text files with common name
filelist = list.files(pattern = "*.txt")

# make a 3 column table of listname,gene,1
res<-lapply(filelist, function(x){
  data.frame(
    set=x,
    geneID=as.character(read.table(x)[,1]),
    val=1)
})

res<-ldply(res)

# turn the 3 column long table into wide
res1<-acast(res,geneID~set,value.var="val",fill=0) 

# force as dataframe
res1<-as.data.frame(res1)

# 1st column must be a name
res1$name=rownames(res1)

# rearrange columns
res2=res1[,c(ncol(res1),1:(ncol(res1)-1))]
head(res2)
colnames(res2) <- c("name", "F_F2IVMvMHco3", "F_MHco18vMHco3","M_F2IVMvF2CTL", "M_F2IVMvMHco3", "M_MHco18vMHco3" )

#make the plot with 5 sets 
upset(res2, mainbar.y.label = "Downregulated Genes", sets = c("F_F2IVMvMHco3", "F_MHco18vMHco3","M_F2IVMvF2CTL", "M_F2IVMvMHco3", "M_MHco18vMHco3"), keep.order = TRUE, nintersects = NA, nsets = 5, text.scale = c(2, 2, 1.5, 1.5, 2, 1.6))


