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
#collmn names are the patient names, and the row names are the gene names
d0 <- DGEList(counts)
d0$counts
d0$samples
table(d0$samples$group)
#pre-processing
d0 <- calcNormFactors(d0)
#Note: calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.
#filtering low expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
snames <- colnames(counts) # Sample names
snames
cultivar <- substr(snames, 1, nchar(snames) - 2) 
time <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
cultivar
time
#our expirements have two factors: cultivar {C, I5, I8} and time {6 and 9}
#Create a new variable “group” that combines cultivar and time

group <- interaction(cultivar, time)
group
plotMDS(d, col = as.numeric(group))
#Voom transformation and calculation of variance weights
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)
#fitting linear models in limma
fit <- lmFit(y, mm)
head(coef(fit))
#Specify which groups to compare:

#Comparison between times 6 and 9 for cultivar I5
contr <- makeContrasts(groupI5.9 - groupI5.6, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
#logFC: log2 fold change of I5.9/I5.6
#AveExpr: Average expression across all samples, in log2 CPM
#t: logFC divided by its standard error
#P.Value: Raw p-value (based on t) from test that logFC differs from 0
#adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
#B: log-odds that gene is DE (arguably less useful than the other columns)
length(which(top.table$adj.P.Val < 0.05))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "time9_v_time6_I5.txt", row.names = F, sep = "\t", quote = F)
