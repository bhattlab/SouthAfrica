# plot metaphlan2 results for ZA data and global data

library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(genefilter)
library(cowplot)

setwd("/Users/Fiona/scg4_fiona/za/za-microbiome/")

pheno_global <- read.table("metadata_final/pheno_global.tsv", sep = "\t", header = T)
za_meta <- read.table("metadata_final/za_meta.tsv", sep = "\t", header = T)

metaphlan_out <- read.table("input_final/merged_abundance_table.txt", sep = "\t", header = T, check.names = F)

# get species
metaphlan_species <- filter(metaphlan_out, grepl("s__", ID) & !grepl("t__", ID))
metaphlan_species$ID <- gsub(".+\\|s__", "", metaphlan_species$ID)
row.names(metaphlan_species) <- metaphlan_species$ID
metaphlan_species$ID <- NULL

# find taxa present in at least 
keep <- data.frame(genefilter(za_species_perc, pOverA(p=0.001, A=0.01)))
keep$tax <- row.names(keep)
keep <- filter(keep, `genefilter.za_species_perc..pOverA.p...0.001..A...0.01..` == T)$tax
za_filt <- za_species_perc[keep, ]

# plot
n_taxa <- 20

# color palette for n taxa
myCols <- colorRampPalette(brewer.pal(12, "Paired"))
my_pal <- myCols(n_taxa)
my_pal <- sample(my_pal)
my_pal[n_taxa + 1] <- "gray"

# find top n taxa
abundance_threshold <- sort(rowSums(metaphlan_species), decreasing = T)[n_taxa]
species_plot <- metaphlan_species[rowSums(metaphlan_species) >= abundance_threshold,]

# add "other" column
species_plot <- rbind(species_plot, t(data.frame("Other" =  100 - colSums(species_plot))))

# melt data frame
species_plot$species <- row.names(species_plot)
species_long <- melt(species_plot, id.vars = "species", variable.name = "sample", value.name = "rel_abundance")

# merge in pheno date
species_meta <- merge(species_long, pheno_global, by = "sample")
# bracken_pheno$label <- labels[match(bracken_pheno$sample, labels$sample_code), "id"]

# set factor level for correct plotting order
# species_meta$species <- factor(species_meta$species, levels = species_meta$species)

# plot in order of decreasing relative abundance of desired taxon
sorted <- filter(species_meta, species == "Prevotella_copri")[, c("sample", "rel_abundance")]
# bracken_pheno$label <- factor(bracken_pheno$label, levels = sorted[rev(order(sorted$rel_abundance)), ]$label)


g <- ggplot(species_meta, aes(x=sample, y=rel_abundance, fill=species)) +
  geom_bar(stat="identity") +
  ylab("Genus-level Relative Abundance") +
  xlab("Sample") +
  scale_fill_manual(values=my_pal) +
  guides(fill = guide_legend(ncol=1, keywidth = 0.125, keyheight = 0.1, default.unit = "inch")) +
  theme_cowplot(12) +
  theme(
    axis.text.x = element_text(size = 4, angle = 90, hjust = 1),
    legend.text = element_text(size = 10)
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))




agin <- plot_bracken(filter(bracken_pheno, site == "Agincourt"), "Agincourt") + theme(legend.position = "none")
soweto <- plot_bracken(filter(bracken_pheno, site == "Soweto"), "Soweto")

za_genus <- plot_grid(agin, soweto, labels = c('A', 'B'), align = "v", ncol = 1, axis = "l", label_size = 12)
ggsave("rplots/supplementary_figures/za_genus.pdf", za_genus, device = "pdf", width = 8, height = 7)
ggsave("rplots/supplementary_figures/za_genus.png", za_genus, device = "png", width = 8, height = 7)