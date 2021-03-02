      
#============================================================================================================================
#
#       This R script will accept the expressio values (such as RPKM) and calculate the coexpression
#
#============================================================================================================================

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one arg ument: if not, return an error
if (length(args)< 1) {
  stop("At least one argument must be supplied (input file).\n USES : Rscript correlation.R tissue.expression ", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.png"
}

  #check if required libraries are installed. Install the required libraries
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!require("gplots")) BiocManager::install("gplots")
  if (!require("gplots")) BiocManager::install("RColorBrewer")


# load the required libraries
library(RColorBrewer) 
library(gplots)

#load the expression data
rpkm <- read.csv(args[1],header=TRUE,row.names=1,check.names = F)
x <- subset(rpkm, select = c(2,3, 4, 5,6)) # change the column number here if there are more less columns are present in dataset
      
#transpose and view datafile
y <- t(x)

      
#define method of correlation analysis: Pearson, Kendall or Spearman
cU<- cor(y, method = c("pearson")) #method can be "pearson", "kendall", "spearman"
write.table(cU, file = "tissue_coexpression_matrix.txt")
rowcol <- as.character(rpkm$color) #choose the column with row color infomration
      
# Prepare the heat map and save as pdf
figname=args[1]
pdf(paste0(args[1], ".pdf"), width = 15, height = 16, family = "Times", pointsize = 16) # defaults to 7 x 7 inches
heatmap.2(cU, key=T,  
keysize = 1, 
key.title = NA, 
trace="none",
col=colorpanel(100, "navy", "white", "firebrick3"), #choose color combination
margin=c(10, 10),
density.info="none",
colRow = rowcol,
colCol = rowcol )
dev.off()
         