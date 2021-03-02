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

1. <bold>COEXPRESSION ANALYSIS: <code> coexpression.R </code> </bold> 
   This R script requires two libraries <code>RColorBrewer</code>  and <code>gplots</code>, and can be run from Ubuntu terminal like this 
   
   <code>Rscript coexpression.R tissue</code>
   
The file <code> tissue </code> contains the expression value iin the form of RPKM (Reads per kilo base per million mapped reads) 

After successful completeion, two files will be generated: This correlation matrix (<code>tissue_coexpression_matrix.txt</code>) and a heatmap (<code>tissue.expression.pdf</code>)

2. <bold> TRANSCRIPTION FACTOR BINDING SITE ANALYSIS: <code> JASPAR.R</code>  </bold> 
  This R script requires following libraries <code>TFBSTools</code>,  <code>JASPAR2020</code>, <code>Biostrings</code>, <code>GenomicRanges</code>, <code>dplyr</code>, <code>ggplot2</code>, and can be run from Ubuntu terminal like this
  
   <code>Rscript JASPAR.R G10Hpro</code>
   
   The file <code> G10Hpro </code> contains the promoter sequence in FASTA format. After successful completeion, a directory  <code>G10Hpro.result </code> with raw data, a table  <code>G10Hpro.result.csv </code> and an image  <code>G10Hpro.png </code> file will be created.  
   
   
   # Citation
   
   Please cite this book chapter if you find these scripts useful for your work
   
   
   
   For any question and comment, please email me -  ssi228@uky.edu

