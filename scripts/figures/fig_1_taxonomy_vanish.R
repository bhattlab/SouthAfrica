library(cowplot)
library(ggpubr)
library(here)
library(tidyverse)

# load data ----
load(here("RData/metadata.RData"))
load(here("RData/za_data.RData"))

# crassphage ----
crass <- za_G[grep("crAss", rownames(za_G)), ]

crass_tot <- crass %>%
  rownames_to_column("crass") %>%
  pivot_longer(cols = -crass, names_to = "sample", values_to = "value") %>%
  left_join(za_meta, by = "sample") %>%
  group_by(sample, site) %>%
  summarise(tot = sum(value)) %>%
  ungroup() %>%
  mutate(present = tot >= 650) %>%
  select(site, present)

fisher.test(table(crass_tot))

crass_tot %>%
  count(site, present)

# panel A: taxonomy barplots ----

sort_by <- "Prevotella"
n_taxa <- 20

# color palette for n taxa
if (!file.exists(here("rds/barplot_pal.rds"))){
  myCols <- colorRampPalette(brewer.pal(9, "Set1"))
  
  barplot_pal <- myCols(n_taxa)
  barplot_pal <- sample(barplot_pal)
  barplot_pal[n_taxa + 1] <- "gray"
  
  saveRDS(barplot_pal, "rds/barplot_pal.rds")
} else {
  barplot_pal <- readRDS(here("rds/barplot_pal.rds"))
}

# find top n taxa
abundance_threshold <- sort(rowSums(za_G_rel), decreasing = T)[n_taxa]
bracken_plot <- za_G_rel[rowSums(za_G_rel) >= abundance_threshold,]

# add "other" column
bracken_plot <- rbind(bracken_plot,
                      t(data.frame("Other" =  1 - colSums(bracken_plot))))

# format for plotting
bracken_plot <- bracken_plot %>%
  rownames_to_column("genus") %>%
  # format labels: remove "miscellaneous", update "environmental Firmicutes"
  mutate(genus = gsub("\\(miscellaneous\\)", "", genus),
         genus = gsub("(Firmicutes): environmental samples",
                      "unclassified \\1", genus))

# pivot to long format
bracken_long <- bracken_plot %>%
  pivot_longer(cols = -genus, names_to = "sample",
               values_to = "rel_abundance") %>%
  # merge pheno data
  left_join(za_pheno[, c("sample", "label_abbrev", "site")], by = "sample") %>%
  group_by(sample) %>%
  # plot in order of decreasing relative abundance of desired taxon
  mutate(top = rel_abundance[genus == sort_by]) %>%
  ungroup() %>%
  mutate(genus = factor(genus, levels = bracken_plot$genus),
         label_abbrev = fct_reorder(label_abbrev, top, .desc = T))

plot_bracken <- function(counts, title){
  g <- ggplot(counts, aes(label_abbrev, rel_abundance * 100, fill = genus)) +
    geom_bar(stat="identity") +
    labs(title = title,
         x = "Participant",
         y = "Relative Abundance (%)",
         fill = "Genus") +
    scale_fill_manual(values = barplot_pal) +
    guides(fill = guide_legend(ncol = 1, keywidth = 0.125, keyheight = 0.1,
                               default.unit = "inch")) +
    theme_cowplot(12) +
    theme(plot.title = element_text(face = "plain", size = 14),
          axis.text.x = element_text(size = 6, angle = 90, hjust = 1,
                                     vjust = 0.5),
          legend.text = element_text(size = 10)) +
    scale_y_continuous(limits = c(0, 100.1), expand = c(0, 0))
  
  return(g)
}

# plot
bbr <- plot_bracken(filter(bracken_long, site == "Bushbuckridge"),
                     "Bushbuckridge") +
  theme(legend.position = "none")

soweto <- plot_bracken(filter(bracken_long, site == "Soweto"), "Soweto")

a <- plot_grid(bbr, soweto, align = "v", ncol = 1, axis = "l")

# panel B: VANISH taxa ----

# pull counts for VANISH taxa
vanish_tax <- data.frame(
  "F" = gsub(".+f__|\\|g__.+", "", vanish_G),
  "G" = gsub(".+g__", "", vanish_G)
)

za_G_VANISH <- za_G_rel[which(rownames(za_G_rel) %in%
                                as.character(vanish_tax$G)), ]
za_F_VANISH <- za_F_rel[which(rownames(za_F_rel) %in%
                                as.character(vanish_tax$`F`)), ]
# long format
counts_long <- za_G_VANISH %>%
  rownames_to_column("G") %>%
  pivot_longer(cols = -G, names_to = "sample", values_to = "relab") %>%
  inner_join(vanish_tax, by = "G") %>%
  inner_join(za_meta, by = "sample") %>%
  mutate(G = fct_reorder(G, relab, .fun = median, .desc = T),
         site = factor(site))

# compare BBR vs SWT means
counts_long %>%
  group_by(G, site) %>%
  summarise(mean = mean(relab)) %>%
  group_by(G) %>%
  summarise("inc_BBR" <- mean[site == "Bushbuckridge"] > mean[site == "Soweto"])

# add pseudo percent
min_perc <- min(counts_long$relab[counts_long$relab > 0])
counts_long$relab[counts_long$relab == 0] <- min_perc / 2

b <- ggplot(counts_long, aes(G, relab * 100, fill = site)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3),
              alpha = 0.75, aes(group = site), size = 1, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  facet_grid(~`F`, scales = "free", space = "free") +
  scale_y_log10(labels = c("0.0001", "0.01", "1", "100"),
                breaks = c(0.0001, 0.01, 1, 100)) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   face = "italic"),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(face = "italic")) +
  scale_fill_manual(values = za_pal) +
  scale_color_manual(values = za_pal) +
  labs(x = "",
       y = "Relative Abundance (%)",
       fill = "") +
  stat_compare_means(label = "p.signif", label.y = 2.5,
                     method = "wilcox.test") +
  background_grid(major = "y")

# figure 1 ----
plot_grid(a, b, labels = c("A", "B"), ncol = 1, rel_heights = c(0.55, 0.45))

ggsave(here("final_plots/figure_1.png"), width = 10, height = 12, bg = "white")
