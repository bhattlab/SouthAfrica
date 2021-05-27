library(plyr)
library(dplyr)
library(data.table)

input <- snakemake@input[[1]]
output <- snakemake@output[[1]]

cazy <- read.table(input, sep = '\t')
# cazy <- read.table("/Users/Fiona/scg4_fiona/prebio2/04_cazy/01_diamond/P1.tsv", sep = '\t')

colnames(cazy) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
cazy_filter <- filter(cazy, evalue < 1e-6)
cazy_filter <- as.data.table(cazy_filter)
cazy_filter_top <- data.frame(cazy_filter[cazy_filter[, .I[which.max(evalue)], by=qseqid]$V1])
cazy_filter_top$sseqid <- gsub('^[^\\|]+\\|', '', cazy_filter_top$sseqid)
s <- strsplit(as.character(cazy_filter_top$sseqid), split = "\\|")
df <- data.frame(sseqid = rep(cazy_filter_top$sseqid, sapply(s, length)), cazyme = unlist(s))

# count table
cazy_counts <- plyr::count(df, "cazyme")
# cazy_counts_wide <- dcast(cazy_counts, ~ cazyme)
cazy_counts$sample <- gsub('.+/|\\.tsv', '', input)
write.table(cazy_counts, output, sep = '\t', row.names = F, quote = F)
