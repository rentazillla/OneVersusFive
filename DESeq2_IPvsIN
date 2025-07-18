# DESeq2 comparison of the IP to IN RNAs

# Input data are treated as the reference level in DESeq2
# setting false discovery rate (alpha) to 0.1
res <- results(dds, contrast = c('Type', 'IP', 'IN'))

# you can call the data to look at it
res

# example of how to perform shrinkage of effect sizes
# improves the accuracy of the padj values
resLFC <- lfcShrink(dds, coef="Type_IP_vs_IN", type = "apeglm")
resLFC

# adding gene symbols from mm.gene_symbols to the resLFC
resLFC <- tibble::rownames_to_column(as.data.frame(resLFC), "GENEID")
resLFC <- left_join(resLFC, mm.gene_symbols, by = "GENEID")

# adding transcripts per kilobase million to the resLFC
length_means <- DataFrame(rowMeans(txi$length))
rpk <- resLFC$baseMean / (length_means$rowMeans.txi.length./1000)
scalingFactor <- sum(rpk)/1e6
resLFC$TPM <- rpk/scalingFactor

# plotting some differential expression data using ggplot2
# results for IP vs IN (all samples)
# we select some importnat microglia genes to highlight how enriched they are
# relative to the bulk genes from the striatum

microglia_genes <-c('Tmem119', 'C1qc', 'Csf1r', 'Selgpl', 
                    'Olfml3', 'Ctss', 'Golm1')


# all genes with log2fc > |1.5| and padj < 0.05
#labels selected microglia genes
myPalette <- colorRampPalette((brewer.pal(11, "Spectral")))
g <- ggplot(resLFC,
            aes(x = log2(TPM), y = log2FoldChange, color = padj,
                   label = SYMBOL)) +
  geom_point() + 
  gghighlight(abs(log2FoldChange) > 1.5 &
                padj < 0.05, 
              unhighlighted_params = list(color = NULL, alpha = 1)) +
  scale_color_gradientn(colors = myPalette(8))

# adding microglia gene annotations
g +
  geom_text_repel(data = subset(resLFC, SYMBOL %in% microglia_genes), 
                                 aes(label = SYMBOL),
                  force = 0.15,
                  nudge_x = 6,
                  direction = 'y',
                  hjust = 0,
                  max.overlaps = Inf)

# using the ggsave function to save plots at any dpi you desire
# ggsave(ggsave("C:/Users/David.SIBCR-1081/Desktop/Figures/IP_vs_IN.png",
#            plot = g,
#             dpi = 300))
