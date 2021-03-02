#================================================================================
#
#       Using "JASPAR2020' packages to detect TFBS
#
#================================================================================

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one arg ument: if not, return an error
if (length(args)< 1) {
  stop("At least one argument must be supplied (input file).\n USES : Rscript JASPAR.R promoter.fasta ", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.png"
}


## Provide the dir name(i.e sub dir) that you want to create under main dir:
output_dir <- file.path('.', file=paste0(args[1], ".result"))

if (!dir.exists(output_dir)){
  dir.create(output_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
} else {
  print("Dir already exists!")
}
sud_dir <- paste0(args[1], ".result")


#check if required libraries are installed. Install the required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("TFBSTools")) BiocManager::install("TFBSTools")
if (!require("JASPAR2020")) BiocManager::install("JASPAR2020")
if (!require("Biostrings")) BiocManager::install("Biostrings")
if (!require("GenomicRanges")) BiocManager::install("GenomicRanges")
if (!require("dplyr")) BiocManager::install("dplyr")
if (!require("ggplot2")) BiocManager::install("ggplot2")


  
## Load Packages
library(TFBSTools)
library(JASPAR2020)
library(Biostrings)
library(GenomicRanges)
library(dplyr)
library(ggplot2)


### For TFBS detection, there are two files needs to be prepared. The first is a PWM (position weighted matrix)
### which has been stored in JASPAR2020, the other is fasta file of the regions you are interested in.
### Make an 'opts' object to help you retrieve the PWM from JASPAR2020
# Query JASPAR database for plants sites:
#getMatrixSet
opts = list()
opts[['tax_group']] <- 'plants'   # other options : "plants", "vertebrates", "insects", "urochor-dat", "nematodes", "fungi"
opts[['matrixtype']] <- 'PWM' # other options : "PFM","PWM","ICM"
PWMatrixList <- getMatrixSet(JASPAR2020, opts)
  

### Make you sequence object
# Load the promoter sequence in fasta format
dnaSet <- readDNAStringSet(args[1], format = 'fasta')
dnaTab <- as.data.frame(dnaSet)


### 'searchSeq()' key function which contains PWM and sequence information. 
#'min.score' setting for restriction. The higher the number the stricter it will be.
#strand = "*" means check the positive and negative strand
j <- 0
for (i in rownames(dnaTab)){
  j <- j+1
  seqName <- i
  seqFas <- as.character(dnaTab[seqName,])
  subject <- DNAString(seqFas)
  siteSetList <- searchSeq(PWMatrixList, subject, seqname = seqName, min.score = '95%', strand = '*')
  if (j == 1){
    write.table(writeGFF3(siteSetList, scoreType = 'relative'), file = 'result_plants_2020_min.score_95.gff3', quote = F,
                sep = '\t', row.names = F)
  } else {
    write.table(writeGFF3(siteSetList, scoreType = 'relative'), file = 'result_plants_2020_min.score_95_not.gff3', quote = F,
                sep = '\t', append = T, col.names = F, row.names = F )
  }
}


#Find scores and remove TFs with no hits
scores <- relScore(siteSetList)
scores <- scores[sapply(scores, length) > 0]
gr_tf <- as(siteSetList[names(scores)], "GRanges")
#Putative binding site with relscore >0.95
gr_tf_filter <- gr_tf[ gr_tf$relScore >0.95 ]

# Summary of overlapping hits.
idx <- unlist(order(mcols(gr_tf_filter)["relScore"], decreasing = TRUE))
as.data.frame(gr_tf_filter[idx]) %>%
  `rownames<-`(NULL) %>%
  select(ID, TF, relScore, absScore, siteSeqs, class) %>% 
  write.csv(.,file = "summary_of_primary_analysis.csv", row.names = F, quote = F) 

# open primary result and custom JASPAR database
df<- read.csv("summary_of_primary_analysis.csv",header=TRUE,check.names = F)
db<- read.csv("database/JASPAR2020",header=TRUE,check.names = F)


#delete duplicate entries
df1 = df[!duplicated(df$ID),]
df2 <- df1[,c(1)]
# frep JASPAR custom file for uniprot IDs
result <- filter(db, grepl(paste(df2, collapse="|"), ID))

#Plot TF binding positions and relative scores
PLOT  <- ggplot(as.data.frame(df1), aes(x = ID, y = relScore, color = class)) +
  geom_label(aes(label = TF)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")
ggsave(PLOT , file=paste0(args[1], ".png"), units="in", family="Arial", width=10, height=5,  dpi = 600, pointsize = 2)

merged <-  merge(df1, result, by.x = 1, by.y = 1, all.x = TRUE)
write.csv(merged,file=paste0(args[1], ".result.csv"), row.names = F, quote = F)



########################################################################################
#
#            Rearrange files
#
#
########################################################################################

if(file.exists("result_plants_2020_min.score_95.gff3") ){
  file.copy("result_plants_2020_min.score_95.gff3",  sud_dir)
} else {
  if(file.exists("result_plants_2020_min.score_95.gff3") ){
    file.copy("result_plants_2020_min.score_95_not.gff3",  sud_dir)
  } else {
    print("Analysis complete. Thanks!")
  }
}

file.copy("summary_of_primary_analysis.csv",  sud_dir)
file.remove("summary_of_primary_analysis.csv")


if(file.exists("result_plants_2020_min.score_95.gff3") ){
  file.remove("result_plants_2020_min.score_95.gff3")
} else {
  if(file.exists("result_plants_2020_min.score_95.gff3") ){
    file.remove("result_plants_2020_min.score_95_not.gff3")
  } else {
    print("Analysis complete. Thanks!")
  }
}


f_tab <-  paste0(args[1], ".result")
f_img <- paste0(args[1], ".png")
message(sprintf("Final result is saved in : %s\n", f_tab))
message(sprintf("Final result is saved in : %s\n", f_img))
message(sprintf("Raw result is saved in : %s\n", sud_dir))
