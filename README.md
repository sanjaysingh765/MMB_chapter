# Methods In Molecular Biology_chapter
In this repository the code for coexpression and promoter analysis is provided. Before running through this tutorial, you must have latest version of R installed.

# R installation
Following command can be use to install and update R on Ubuntu

#Update and Install R
<code>sudo apt-get update
sudo apt-get install r-base r-base-dev
sudo apt-get upgrade</code>

#Know R version
<code>R --version</code>



This script provides an overview of commands for the analysis pipeline of protein complexes from SEC-Evosep-SWATH data. For quantitative data extraction the OpenSWATH Workflow (OSW) (Rost et al. 2014) is applied and subsequent statiscal scoring of co-elution proteins with the R-package CCprofiler. To conduct the trial experiment, download the *.mzXML data, the spectral assay library, and iRT file from ProteomeXchange (data set identifier: XXXX). Installation instruction for applied software tools and detailed instructions about commands and parameters are documented on http://openswath.org/en/latest/ .

