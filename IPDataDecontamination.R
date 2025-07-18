# here we correct the IP gene quantification based on possible contamination arising from 
# any IN material that may have entered the IP samples prioir to sending them off for sequencing
# see our previous paper for the reasoning behind this argument as well as a separate example
# using some simulated RNASeq data
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9015218/#appsec1

# get gene counts for each IP sample (one mouse = one sample)
counts <- as.data.frame(txi$counts) %>%
  dplyr::select(starts_with(c("IP"))) %>%
  rownames_to_column("GENEID")


# get the input values for the genes to correct for contamination
inputs <- as.data.frame(txi$counts) %>%
  dplyr::select(starts_with(c("IN"))) %>%
  rownames_to_column("GENEID")


# sort counts and inputs by sample name
names.IN <-str_sort(names[grep('IN', names)])
names.IP <- str_sort(names[grep('IP', names)])
inputs <- inputs[, c("GENEID", names.IN)]
counts <- counts[, c("GENEID", names.IP)]

# calculate percent contamination of each gene in IP data
# from counts in the input data

# save counts as the cleaned IP data
write.csv(counts, 'IP_Counts_Minus_IN.csv', na = 'NA')

# we don't have all paired IP and IN samples so we will use the select paired
# samples for this calculation

paired <- c("4359", "4363", "4694", "4701", "4704", "4837", "4932", "4933", 
            "4934", "4936", "4937", "4939", "4941", "4973", "4976", "4993",
            "4999", "5004", "5157", "5189", "5568", "5609")

# below is the experimental RNA yield in ng from each IP sample above
ribotag_rna <- c(196.68, 120.56, 125.4, 126.28, 154.88, 127.6, 163.68, 96.8,
                 162.8, 118.8, 92.4, 1.18, 1.29, 154, 105.6, 1.32, 1.50,
                 119.68, 242)

# negative controls below were treated as ribotag animals in extracting IP
# rna but did not express HA tagged ribosomes leading to a very low
# RNA output; these represent average contamination of samples from the IN
# portion of the RNA in ng
negative_controls <- c(1.21, 1.19)

# approximately 20% of the input RNA can be argued to be contaminating
# each IP sample

est_contamination <- mean(mean(negative_controls)/ribotag_rna)

# determine counts of each gene from IN to subtract from IP
# Female Fent Five
for (col in counts[, 3:6]) {
  col = col - (est_contamination * rowMeans(inputs[, 3:5]))
}

# Female Fent One
for (col in counts[, 7:11]) {
  col = col - (est_contamination * rowMeans(inputs[, 6:8]))
}

# Female Sal Five
for (col in counts[, 12:13]) {
  col = col - (est_contamination * rowMeans(inputs[, 9:13]))
}

#Female Sal One
for (col in counts[, 14:18]) {
  col = col - (est_contamination * rowMeans(inputs[, 14:16]))
}

# Male Fent Five
for (col in counts[, 20:22]) {
  col = col - (est_contamination * rowMeans(inputs[, 18:19]))
}

# Male Fent One
for (col in counts[, 23:26]) {
  col = col - (est_contamination * rowMeans(inputs[, 20:23]))
}

# Male Sal Five
for (col in counts[, 27:28]) {
  col = col - (est_contamination * rowMeans(inputs[, 24:27]))
}

# Male Sal One
for (col in counts[, 29:33]) {
  col = col - (est_contamination * rowMeans(inputs[, 28:32]))
}

# make all counts that are negative equal to 0
# and make all counts integers
for (i in 2:33) {
  tmp <- counts[, i]
  for (j in 1:length(tmp)) {
    if (tmp[j] < 0) {
      tmp[j] <- 0 # set negative counts to 0
    }
  }
  counts[, i] <- round(tmp) # counts must be integers
}
