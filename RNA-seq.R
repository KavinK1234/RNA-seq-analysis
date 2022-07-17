#in this project, we will analyze RNA sequencing data from mice breast tissues demonstrating the use of edgeR-limma package
#we will import, organize, filter, and normalize the data to assess differential expression and derive a set of biomarkers
#Sheridan et al. (2015) (Sheridan et al. 2015) and consists of three cell populations (basal, luminal progenitor (LP) and mature luminal (ML)) sorted from the mammary glands of female virgin mice, each profiled in triplicate.


#step 7 though 8 only once
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Glimma")
install.packages("R.utils")
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
#downloading the counts data from github
counts <- read.delim("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/thursday/all_counts.txt")
head(counts)
d0 <- DGEList(counts)
d0$counts
d0$samples
