library(cowplot)
library(ggpubr)
library(here)
library(RColorBrewer)
library(reshape2)
library(tidyverse)

# load data ----
load(here("RData/metadata.RData"))

# nanopore figure ----
labels <- read.table(here("input_final/za_labels.tsv"), sep = "\t", header = T)

# nanopore mags
nmags <- read.table(here("nanopore/input_final/nanopore_gtdbtk.bac120.summary.tsv"),
                            sep = "\t", header = T)

# short read mag data
gtdbtk_bac <- read.table(here("input_final/mags/gtdbtk.bac120.summary.tsv"),
                         sep = "\t", header = T)
gtdbtk_arc <- read.table(here("input_final/mags/gtdbtk.ar122.summary.tsv"),
                         sep = "\t", header = T)
gtdbtk <- rbind(gtdbtk_bac, gtdbtk_arc)
gtdbtk <- gtdbtk %>%
  separate(user_genome, into = c("Sample", "Bin"), sep = "_", extra = "merge")

# binning dastool pipeline output
mags <- read.table(here("input_final/mags/binning_table_all_simple.tsv"),
                   sep = "\t", header = T, comment.char = "")

# format mag info
mags <- mags %>%
  filter(Bin != "unbinned") %>%
  mutate(Quality = ifelse(Completeness > 90 & Contamination < 5 & tRNA >= 18 &
                            rna.16S > 0 & rna.23S > 0 & rna.5S > 0,
                          "High-quality",
                          ifelse(Completeness >= 50 & Contamination <10,
                                 "Medium-quality", "Low-quality")),
         Quality = factor(Quality, levels = c("Low-quality", "Medium-quality",
                                              "High-quality"))) %>%
  left_join(za_meta, by = c("Sample" = "sample")) %>%
  left_join(gtdbtk[, c("Sample", "Bin", "classification")],
            by = c("Sample", "Bin"))

mags_mqplus <- mags %>%
  filter(Quality != "Low-quality")

# gtdb classifications ----
for (lvl in c("S", "G")){
  bracken <- read.table(here(paste0("input_final/kraken/gtdb_r95/bracken_", lvl, ".txt")),
                        sep = "\t", header = T)
  rownames(bracken) <- gsub("[a-z]__", "", rownames(bracken))
  
  # keep relevant samples
  bracken <- bracken[, za_pheno$sample]
  
  # rm zero features
  bracken <- bracken[rowSums(bracken) > 0, ]
  
  bracken_rel <- sweep(bracken, 2, colSums(bracken), FUN = "/")
  
  assign(paste0("bracken_", lvl), bracken)
  assign(paste0("bracken_", lvl, "_rel"), bracken_rel)
}

plot_nanopore <- function(count_data, sample, n_taxa, rank_label){
  
  # color palette for n taxa
  myCols <- colorRampPalette(brewer.pal(12, "Paired"))
  global_pal <- myCols(n_taxa)
  global_pal <- sample(global_pal)
  global_pal[n_taxa + 1] <- "gray"
  
  # find top n taxa
  bracken_plot <- count_data %>%
    select(one_of(sample)) %>%
    rownames_to_column("taxon") %>%
    pivot_longer(cols = -taxon, names_to = "sample",
                 values_to = "rel_abundance") %>%
    arrange(desc(rel_abundance)) %>%
    head(n_taxa)
  
  # add "other" column
  bracken_plot <- rbind(bracken_plot,
                        data.frame(taxon = "Other", sample = sample,
                                   rel_abundance =  1 - sum(bracken_plot$rel_abundance)))
  
  # merge in pheno date
  bracken_plot$label <- labels[match(bracken_plot$sample, labels$sample), "label"]

  # designate whether taxon was assembled into a MAG or nMAG
  pattern <- ifelse(rank_label == "Genus", ".+g__|;s__.*", ".+s__")
  
  nmags_filt <- nmags %>%
    mutate(Sample = gsub("^n|_.+", "", user_genome),
           taxon = gsub(pattern, "", classification)) %>%
    select(Sample, classification, taxon) %>%
    filter(Sample == !!sample)
  
  mags_filt <- mags_mqplus %>%
    mutate(taxon = gsub(pattern, "", classification)) %>%
    select(Sample, taxon) %>%
    filter(Sample == sample)
  
  bracken_plot <- bracken_plot %>%
    mutate(taxon_label = ifelse(taxon %in% mags_filt$taxon,
                                paste0(taxon, "*"),
                                taxon),
           taxon_label = ifelse(taxon %in% nmags_filt$taxon,
                                paste0(taxon_label, "†"),
                                taxon_label),
           taxon = taxon_label,
           taxon = factor(taxon, levels = unique(taxon)))
  
  ggplot(bracken_plot, aes(label, rel_abundance * 100, fill = taxon)) +
    geom_bar(stat="identity") +
    labs(
      x = "",
      y = paste(rank_label, "relative abundance (%)"),
      fill = rank_label
    ) +
    scale_fill_manual(values = global_pal) +
    guides(fill = guide_legend(ncol=1, keywidth = 0.125, keyheight = 0.1,
                               default.unit = "inch")) +
    theme_cowplot() +
    theme(
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0))
}

nmag_samples <- c("C27", "C29", "C33")
nmag_species <- lapply(nmag_samples, function(x){plot_nanopore(bracken_S_rel,
                                                               x,30, "Species")})
nmag_genus <- lapply(nmag_samples, function(x){plot_nanopore(bracken_G_rel,
                                                             x, 30, "Genus")})

a1 <- nmag_species[[1]] + ggtitle("Bushbuckridge 106")
a2 <- nmag_genus[[1]]

b1 <- nmag_species[[2]] + ggtitle("Bushbuckridge 108")
b2 <- nmag_genus[[2]]

c1 <- nmag_species[[3]] + ggtitle("Bushbuckridge 113")
c2 <- nmag_genus[[3]]

plot_grid(a1, a2, b1, b2, c1, c2,
          ncol = 2, align = "hv", axis = "bt",
          labels = c("A", "", "B", "", "C", ""))

ggsave(here("final_plots/supplementary/figure_S16_za_nanopore_taxa2.png"),
       width = 8.5, height = 15)

# # check labels
# s <- "C33"
# nanopore_bins %>%
#   filter(Sample == s) %>%
#   arrange(Final.Class) %>%
#   pull(Final.Class) %>%
#   unique()
# 
# illumina_bins %>%
#   filter(Sample == s) %>%
#   arrange(Final.Class)
# 
