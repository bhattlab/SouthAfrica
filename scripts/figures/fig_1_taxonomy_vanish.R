library(cowplot)
library(ggpubr)
library(here)
library(reshape2)
library(tidyverse)

## load data ----
source(here("scripts/load_data.R"))

labs <- read.table(here("input_final/za_labels.tsv"), sep = "\t", header = T)

## Figure 1A: taxonomy barplots ----
sort_by <- "Prevotella"
n_taxa <- 20

# color palette for n taxa
if (!file.exists(here("rds/global_pal.rds"))){
  myCols <- colorRampPalette(brewer.pal(9, "Set1"))
  barplot_pal <- myCols(n_taxa)
  barplot_pal <- sample(barplot_pal)
  barplot_pal[n_taxa + 1] <- "gray"
  saveRDS(barplot_pal, "rds/global_pal.rds")
} else {
  barplot_pal <- readRDS("rds/global_pal.rds")
}

# find top n taxa
abundance_threshold <- sort(rowSums(za_G_rel), decreasing = T)[n_taxa]
bracken_plot <- za_G_rel[rowSums(za_G_rel) >= abundance_threshold,]

# add "other" column
bracken_plot <- rbind(bracken_plot, t(data.frame("Other" =  1 - colSums(bracken_plot))))

# melt data frame
bracken_plot$Genus <- row.names(bracken_plot)

# remove annoying "miscellaneous" labels that take up space
bracken_plot$Genus <- gsub("\\(miscellaneous\\)", "", bracken_plot$Genus)
bracken_long <- melt(bracken_plot, id.vars = "Genus", variable.name = "sample", value.name = "rel_abundance")

# merge in pheno date
bracken_pheno <- merge(bracken_long, za_meta, by = "sample")
bracken_pheno$label <- labs[match(bracken_pheno$sample, labs$sample), "label"]

# set factor level for correct plotting order
bracken_pheno$Genus <- factor(bracken_pheno$Genus, levels = bracken_plot$Genus)

# shorten labels
bracken_pheno$label <- gsub("bushbuckridge_", "BBR", bracken_pheno$label)
bracken_pheno$label <- gsub("soweto_", "SWT", bracken_pheno$label)

# plot in order of decreasing relative abundance of desired taxon
sorted <- bracken_pheno %>% filter(Genus == sort_by) %>% arrange(desc(rel_abundance)) %>% pull(label)
bracken_pheno$label <- factor(bracken_pheno$label, levels = sorted)

plot_bracken <- function(counts, title){
  g <- ggplot(counts, aes(x=label, y=rel_abundance * 100, fill=Genus)) +
    geom_bar(stat="identity") +
    labs(
      title = title,
      x = "Participant",
      y = "Relative Abundance (%)"
    ) +
    scale_fill_manual(values = barplot_pal) +
    guides(fill = guide_legend(ncol=1, keywidth = 0.125, keyheight = 0.1, default.unit = "inch")) +
    theme_cowplot(12) +
    theme(
      plot.title = element_text(face = "plain", size = 14),
      axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
      legend.text = element_text(size = 10)
    ) +
    scale_y_continuous(limits = c(0, 100.1), expand = c(0, 0))
  
  return(g)
}

agin <- plot_bracken(filter(bracken_pheno, site == "Bushbuckridge"), "Bushbuckridge") + theme(legend.position = "none")
soweto <- plot_bracken(filter(bracken_pheno, site == "Soweto"), "Soweto")

a <- plot_grid(agin, soweto, align = "v", ncol = 1, axis = "l")


## Figure 1B: VANISH taxa ----
vanish_tax <- data.frame(
  "F" = gsub(".+f__|\\|g__.+", "", vanish_G),
  "G" = gsub(".+g__", "", vanish_G)
)

# za_G_VANISH <- za_G_pseudo_rel[which(rownames(za_G_pseudo_rel) %in% as.character(vanish_tax$G)), ]
# za_F_VANISH <- za_F_pseudo_rel[which(rownames(za_F_pseudo_rel) %in% as.character(vanish_tax$`F`)), ]

za_G_VANISH <- za_G_rel[which(rownames(za_G_rel) %in% as.character(vanish_tax$G)), ]
za_F_VANISH <- za_F_rel[which(rownames(za_F_rel) %in% as.character(vanish_tax$`F`)), ]

counts_long <- reshape2::melt(za_G_VANISH %>% rownames_to_column(var = "G"),
                              id.vars = "G",
                              variable.name = "sample",
                              value.name = "relab")

counts_long$G <- factor(counts_long$G, levels = rev(unique(counts_long$G)))
counts_long <- merge(counts_long, za_meta, by.x = "sample", by.y = "sample")

# facet by family
counts_long$F <- vanish_tax[match(counts_long$G, vanish_tax$G), "F"]

# compare BBR vs SWT means
counts_long %>%
  group_by(G, site) %>%
  summarise(mean = mean(relab)) %>%
  group_by(G) %>%
  summarise("inc_BBR" <- mean[site == "Bushbuckridge"] > mean[site == "Soweto"])

# plot in decreasing order
levels <- counts_long %>%
  group_by(G) %>%
  summarise(mean_relab = mean(relab)) %>%
  arrange(-mean_relab) %>%
  pull(G)
counts_long$G <- factor(counts_long$G, levels = levels)

counts_long$relab[counts_long$relab == 0] <- min(counts_long$relab[counts_long$relab > 0]) / 2

b <- ggplot(counts_long, aes(G, relab * 100, fill = site)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.75, aes(group = site), size = 1, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  facet_grid(~ F, scales = "free", space = "free") +
  scale_y_log10(labels = c("0.0001", "0.01", "1", "100"), breaks = c(0.0001, 0.01, 1, 100)) +
  theme_cowplot() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(face = "italic"),
  ) +
  scale_fill_manual(values = za_pal) +
  scale_color_manual(values = za_pal) +
  labs(
    x = "",
    y = "Relative Abundance (%)",
    fill = ""
  ) +
  stat_compare_means(label = "p.signif", label.y = 2.5) +
  background_grid(major = "y")


## Figure 1 ----
plot_grid(a, b, labels = c("A", "B"), ncol = 1, rel_heights = c(0.55, 0.45))
ggsave(here("final_plots/figure_1.png"), width = 10, height = 12)
