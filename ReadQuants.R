# reading quantified RNAseq data into R

library(tximport)
library(readr)
library(tidyverse)
library(ensembldb)
library(RMariaDB)
library(AnnotationHub)
library(DESeq2)
library(apeglm)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(gplots)
library(RColorBrewer)
library(gghighlight)
library(EDASeq)
library(DEGreport)
library(RNAseqQC)
library(clusterProfiler)
library(thematic)

# create table for the samples of the RNAseq
names <- list.files('Fastqs/', pattern = '._R1\\.fastq.gz', full.names = FALSE)
names <- str_remove(names, '.fastq.gz')

# create a data frame for all the information shown in the sample file names
# this is the sample information data which will be combined with the
# count data on each gene for each sample
samples <- data.frame(do.call(rbind, strsplit(names, "_")))
colnames(samples) <- c('Type', 'Sex', 'Drug', 'Cycle',
                       'Mouse', 'SampleNum', 'SampleName')
samples <- dplyr::select(samples, -7)
samples$SampleName <- names

# save samples as csv
write.csv(samples, file = 'SampleInformation.csv')

# getting all the quantification files from Salmon
# file paths
files <- file.path('quants', str_c(names, "_quant"), 'quant.sf')
# set names
names(files) <- paste0(samples$SampleName)

# get the gene annotations for gene IDs based on the Mouse reference
txdb <- makeTxDbFromEnsembl(organism = "Mus musculus", release = 111)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- ensembldb::select(txdb, k, "GENEID", "TXNAME")

# map the counts from the salmon data based on the mouse reference
# and load the data into r
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)
names(txi)
head(txi$counts)

# saving the count data for each sample
# can be combined with sample information table to rebuild dds
write.csv(txi$counts, file = 'CountData.csv')
write.csv(txi$length, file = 'LengthData.csv')
