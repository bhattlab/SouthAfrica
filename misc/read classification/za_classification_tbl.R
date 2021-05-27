library(dplyr)
library(ggplot2)
library(cowplot)

readcounts <- read.table("/Users/tamburif/za-microbiome/input_final/readcounts/readcounts_preproc_za.tsv", sep = "\t", header = T)

read_class <- read.table("/Users/tamburif/Desktop/classified_all.txt", sep = "\t", header = F, col.names = c("sample", "rank", "reads"))

# multiply by two
# read_class$reads <- read_class$reads * 2

read_class <- read_class %>%
  mutate(
    reads_total = readcounts[match(sample, readcounts$Sample), "host_removed_reads"],
    reads = ifelse(rank == "U", reads_total - reads, reads),
    rank = ifelse(rank == "U", "D", rank),
    perc_classified = reads / reads_total,
    perc_unclassified = (reads_total - reads) / reads_total
    )

read_class$rank <- factor(read_class$rank, levels = c("D", "P", "C", "O", "F", "G", "S"))

ggplot(read_class, aes(rank, perc_classified, fill = rank)) + 
  geom_boxplot() +
  theme_cowplot()

tbl <- read_class %>%
  group_by(rank) %>%
  summarise(
    mean_perc_classified = round(mean(perc_classified), 3),
    range_perc_classified = paste(round(range(perc_classified), 3), collapse = "-"),
    mean_perc_unclassified = round(mean(perc_unclassified), 3),
    range_perc_unclassified = paste(round(range(perc_unclassified), 3), collapse = "-"),
    )
names(tbl) <- c("Taxonomic rank", "Mean classified reads (%)", "Range classified reads (%)", "Mean unclassified reads (%)", "Range unclassified reads (%)")

write.table(tbl, "/Users/tamburif/Desktop/2020-9-18_read_classification_perc.txt", sep = "\t", quote = F, row.names = F)
