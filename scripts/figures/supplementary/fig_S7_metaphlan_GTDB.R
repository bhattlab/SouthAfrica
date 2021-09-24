library(cowplot)
library(ggpubr)
library(here)
library(RColorBrewer)
library(reshape2)
library(scales)
library(tidyverse)
library(vegan)

# load data ----
load(here("RData/metadata.RData"))
load(here("RData/za_data.RData"))

# read gtdb data ---
for (lvl in c("S", "G", "P")){
  bracken <- read.table(here(paste0("input_final/kraken/gtdb_r95/bracken_", lvl, ".txt")),
                        sep = "\t", header = T)
  rownames(bracken) <- gsub("s__", "", rownames(bracken))
  
  # keep relevant samples
  bracken <- bracken[, za_pheno$sample]
  
  # rm zero features
  bracken <- bracken[rowSums(bracken) > 0, ]
  
  bracken_rel <- sweep(bracken, 2, colSums(bracken), FUN = "/")
  
  assign(paste0("bracken_", lvl), bracken)
  assign(paste0("bracken_", lvl, "_rel"), bracken_rel)
}

# function ----

plot_tax <- function(counts, ntax = 20, legend_ncol = 4){
  cts_long <- counts %>%
    rownames_to_column("feature") %>%
    pivot_longer(cols = !(feature), names_to = "sample", values_to = "relab") %>%
    group_by(feature) %>%
    mutate(mean_relab = mean(relab)) %>%
    arrange(-mean_relab) %>%
    ungroup() %>%
    filter(dense_rank(dplyr::desc(mean_relab)) %in% 1:ntax) %>%
    select(-mean_relab)
  
  other <- cts_long %>%
    group_by(sample) %>%
    summarise(relab = 100 - sum(relab)) %>%
    mutate(feature = "Other")
  
  cts_long <- rbind(cts_long, other) %>%
    left_join(za_pheno, by = "sample")
  
  # plot features by relab
  lvls_y <- cts_long %>%
    group_by(feature) %>%
    summarise(mean = mean(relab)) %>%
    arrange(-mean) %>%
    filter(feature != "Other") %>%
    pull(feature) %>%
    as.character()
  
  cts_long$feature <- factor(cts_long$feature, levels = c(lvls_y, "Other"))
  
  # plot samples by top feature
  lvls_x <- cts_long %>%
    filter(feature == lvls_y[1]) %>%
    arrange(-relab) %>%
    pull(sample)
  
  cts_long$sample <- factor(cts_long$sample, levels = lvls_x)
  
  # palette
  myCols <- colorRampPalette(brewer.pal(12, "Paired"))
  my_pal <- myCols(ntax)
  my_pal <- sample(my_pal)
  my_pal[ntax + 1] <- "#BBBBBB"
  
  # plot
  ggplot(cts_long, aes(sample, relab, fill = feature)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = my_pal) +
    labs(x = "Sample",
         y = "Relative abundance (%)",
         fill = "") +
    guides(fill = guide_legend(ncol = legend_ncol)) +
    facet_grid(~site, scales = "free_x", space = "free") +
    theme_cowplot(12) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom")
}

plot_richness <- function(counts, cohort){
  shan <- diversity(counts, index = "shannon", MARGIN = 2)
  
  slong <- melt(shan) %>%
    rownames_to_column("sample") %>%
    mutate(cohort = !!cohort)
}

process_metaphlan <- function(m_out, rank, rank_filter_stop){
  
  rank_prefix <- paste0(rank, "__")
  
  m_out_filt <- m_out %>%
    filter(grepl(rank_prefix, clade_name), !grepl(rank_filter_stop, clade_name)) %>%
    dplyr::select(-NCBI_tax_id) %>%
    column_to_rownames("clade_name")
  
  # keep relevant samples
  m_out_filt <- m_out_filt[, za_meta$sample]
  
  # rm zero features
  m_out_filt <- m_out_filt[rowSums(m_out_filt) > 0, ]
  
  # format feature names
  rownames(m_out_filt) <- make.names(gsub(".+g__", "", rownames(m_out_filt)),
                                     unique = T)
  
  return(m_out_filt)
}

# metaphlan plot ----
m_out <- read.table(here("input_final/metaphlan/metaphlan_merge.txt"),
                    sep = "\t", header = T)

colnames(m_out) <- gsub("_concat.+", "", colnames(m_out))

m_out_g <- process_metaphlan(m_out, "g", "s__")
m_out_s <- process_metaphlan(m_out, "s", "placeholder")

plot_m <- m_out_s
rownames(plot_m) <- make.names(gsub(".+s__", "", rownames(plot_m)), unique = T)

a <- plot_tax(plot_m)

# gtdb plot ----
cts <- bracken_G_rel
rownames(cts) <- gsub("g__", "", rownames(cts))

b <- plot_tax(cts * 100, legend_ncol = 6)

# genbank plot ----
# c <- plot_tax(za_S_rel * 100)

# plot shannon diversity ----

m <- plot_richness(m_out_s/100, "Metaphlan")
gb <- plot_richness(za_S_rel, "Genbank")
gtdb <- plot_richness(bracken_S_rel, "GTDB")

all <- rbind(m, gb, gtdb)

d <- ggplot(all, aes(cohort, value)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox", label = "p.signif", 
                     comparisons = list(c("Genbank", "GTDB"), c("Metaphlan", "GTDB")),
                     vjust = 0.5) +
  theme_cowplot() +
  background_grid(major = "y") +
  labs(x = "Database",
       y = "Shannon diversity")

plot_grid(a, b, plot_grid(d, NULL, nrow = 1),
          ncol = 1, labels = c("A", "B", "C"),
          rel_heights = c(0.36, 0.34, 0.3))

ggsave(here("final_plots/supplementary/figure_S7_metaphlan_gtdb.png"),
       width = 8.5, height = 12)

# stats tests
# Genbank vs GTDB
wilcox.test(value ~ cohort, rbind(gb, gtdb))

# Metaphlan vs GTDB
wilcox.test(value ~ cohort, rbind(m, gtdb))

