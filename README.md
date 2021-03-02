# Methods In Molecular Biology_chapter
In this repository the code for coexpression and promoter analysis is provided. Before running through this tutorial, you must have latest version of R installed.

# R installation
Following command can be use to install and update R on Ubuntu

#Update and Install R

<code>sudo apt-get update</code>

<code>sudo apt-get install r-base r-base-dev</code>

<code>sudo apt-get upgrade</code>

#Know R version

<code>R --version</code>

# Running the Scripts

1. <bold>coexpression analysis: coexpression.R </bold> 
   This R script requires two libraries <code>RColorBrewer</code>  and <code>gplots</code>, can be run from Ubuntu terminal like this 
   
   <code>Rscript coexpression.R tissue</code>

After successful completeion, two files will be generated: This correlation matrix (<code>tissue_coexpression_matrix.txt</code>) and a heatmap (<code>tissue.expression.pdf</code>)

2. <bold> Transcription factor binding sites  analysis: JASPAR.R </bold> 
  This R script requires two libraries <code>TFBSTools</code>,  <code>JASPAR2020</code>, <code>Biostrings</code>, <code>GenomicRanges</code>, <code>dplyr</code>, <code>ggplot2</code>

