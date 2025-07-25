# OneVersusFive
A collection of R and Python notebooks used for analyzing RNAseq and mouse behavioral data

## General steps in the RNAseq quantification process

1. [Quantify](https://github.com/rentazillla/OneVersusFive/blob/main/ReadQuants.R) reads using Salmon (https://combine-lab.github.io/salmon/). This will produce quantification files that can then be read into R for data analysis. Note that using Salmon requires making decisions about which transcriptome to use for a basis comparison to your data, as well as other considerations. For this data we compared the RNAseq reads to the mouse transcriptome available here: https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus_c57bl6nj/cdna/

The Salmon website gives extensive documentation on usage and behavior.

2. [Read](https://github.com/rentazillla/OneVersusFive/blob/main/ReadQuants.R) the quantified RNASeq data from Salmon into R.

3. [Perform](https://github.com/rentazillla/OneVersusFive/blob/main/QCTesting.R) some basic quality control checks to make sure your quantification was successful and your data is ready for analysis. See https://cran.r-project.org/web/packages/RNAseqQC/vignettes/introduction.html for the examples presented. Note that in this data set there are two sets of data for the microglia specific genes (the immunoprecipitated RNAs or IP) and the rest of the bulk RNAs from the striatum that were bulk sequenced (the input or IN). Across this work we generally focus on the IP genes and compare their expression to that in the IN genes, as we were interested in microglial gene changes between one and five cycles of opioid withdrawal.

4. [Apply](https://github.com/rentazillla/OneVersusFive/blob/main/DESeq2_IPvsIN.R) differential expression analysis using the DESeq2 package. We compare IP to IN and determine which genes are differentially expressed in microglia relative to the bulk striatum RNAs.

5. [Clean](https://github.com/rentazillla/OneVersusFive/blob/main/IPDataDecontamination.R) the IP data by using an average calculation of IN RNA potentially present in IP samples before sending off for NGS.

6. [Prepare](https://github.com/rentazillla/OneVersusFive/blob/main/IPGroupComparisons) differential expression comparisons between the various experimental groups. Note that Fent = animals treated with Fentanyl, Sal = animals treated with saline, Five = total of five cycles of withdrawal, One = one experience of withdrawal.
Fent Five vs Fent One
Sal Five vs Sal One
Fent Five vs Sal Five
Fent One vs Sal One
