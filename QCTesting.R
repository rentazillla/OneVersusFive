# perform quality contol testing following reading in the data to R

library("RNAseqQC")
library("DESeq2")
library("ensembldb")
library("dplyr")
library("ggplot2")
library("purrr")
library("tidyr")
library("tibble")
library("magrittr")

# create DESeqDataSet to input into differential expression analysis
# using DESeq2
dds <- DESeqDataSetFromTximport(txi, samples, ~Type)

# adding the gene names to the dds data
# mouse ensembl annotation for release 111
ah <- AnnotationHub()
endb.mmusculus <- ah[["AH116340"]]

# getting gene symbols from geneids
mm.genes <- tx2gene$GENEID
mm.gene_symbols <- ensembldb::select(endb.mmusculus, 
    keys= tx2gene$GENEID,
    keytype = "GENEID", columns = c("SYMBOL","GENEID"))

# setting factors of the highest-level comparison treatment
# Input vs IP samples, Input vs Undetermined samples
dds$Type <- relevel(dds$Type, ref = 'IN')

# running DESeq2
dds <- DESeq(dds)

# examples of some QC parameters you may want to check

# plotting total quants of the data
plot_total_counts(dds)

# plotting library complexity 
plot_library_complexity(dds)

# plotting the types of genes found by their biological function
plot_biotypes(dds)
