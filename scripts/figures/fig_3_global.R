library(cowplot)
library(DESeq2)
library(ggrepel)
library(ggpubr)
library(here)
library(reshape2)
library(scales)
library(tidyverse)
library(vegan)
library(metagenomeSeq)

# load data ----
source(here("scripts/load_data.R"))

# rarefy ----
# rare_lvl <- 1e6
# global_S_rare <- rrarefy(t(global_S), rare_lvl)
# global_S_rare <- data.frame(t(global_S_rare))
# 
# mr <- newMRexperiment(global_S_rare)
# p <- cumNormStatFast(mr)
# mr_css <- cumNorm(mr, p = p)
# global_S_rare_css <- MRcounts(mr_css, norm = T, log = T)
# global_S_rare_css <- data.frame(global_S_rare_css)
# names(global_S_rare_css) <- gsub("^X", "", names(global_S_rare_css))
# 
# global_S_rare_rel <- sweep(global_S_rare, 2, colSums(global_S_rare), FUN = "/")

# panel A: cohort summary plot ----
cohort_size <- pheno_global %>%
  group_by(site2) %>%
  tally()

cohort_plot <- ggplot(cohort_size, aes(site2, n, fill = site2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = global_pal) +
  scale_y_continuous(breaks = seq(0, 150, 25)) +
  theme_cowplot(12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = "",
    y = "Participants"
  ) +
  background_grid(major = "y")

# panel B: MDS ----

# bray-curtis distance
# vare_dis <- vegdist(t(global_S_rare_css), method = "bray")
vare_dis <- vegdist(t(global_S_css), method = "bray")

# permanova
meta <- pheno_global %>%
  filter(sample %in% names(global_S_rel))

meta <- meta[match(names(global_S_rel), meta$sample), ]

adonis2(vare_dis ~ site2, data = meta)

dispersion <- betadisper(vare_dis, group = meta$site2)
permutest(dispersion)

# calculate mds
mds <- cmdscale(vare_dis, eig = TRUE, x.ret = TRUE)

# calculate weighted species score
mds_values <- mds$points

# calculate percentage of variation per mds axis
mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)

# mds scores
mds_data <- data.frame(sample = rownames(mds_values),
                       x = mds_values[,1],
                       y = mds_values[,2])

# merge pheno data
mds_meta <- merge(mds_data, pheno_global, by = "sample")

# flip axes for consistency
mds_meta$x <- mds_meta$x * -1
mds_meta$y <- mds_meta$y * -1

global_mds <- ggplot(mds_meta, aes(x, y, color = site2)) +
  geom_point(size = 2, alpha = 0.75) +
  scale_color_manual(values = global_pal) +
  labs(x = paste("MDS 1 (", mds_var_per[1], "%)", sep = ""),
       y = paste("MDS 2 (", mds_var_per[2], "%)", sep = ""),
       color = " ") +
  theme_cowplot() +
  theme(legend.position = "top",
        legend.justification = "center") +
  coord_fixed() +
  background_grid() +
  guides(color = guide_legend(nrow = 1))

global_mds_ci <- global_mds +
  stat_ellipse(aes(color = site2), type = 't', size = 0.5, alpha = 0.75,
               show.legend = F)

# taxonomy plots (similar to Smits 2017) ----

# identify taxa correlated with mds 1 and 2
global_F_rel_filt <- global_F_rel[genefilter(global_F_rel, pOverA(0.1, 0.01)), ]

global_F_t <- t(global_F_rel_filt)
mds1_values <- merge(mds_meta[, c("sample", "x", "y")], global_F_t,
                     by.x = "sample", by.y = "row.names", all = T) %>%
  column_to_rownames(var = "sample")

mds1_cor <- cor(mds1_values, method = "spearman")
mds1_cor_long <- melt(mds1_cor, value.name = "cor")

n_features <- 4

features_x <- mds1_cor_long %>%
  filter(Var1 == "x" & Var2 != "x") %>%
  arrange(-abs(cor)) %>%
  head(n_features) %>%
  pull(Var2) %>%
  as.character()

# figure legend rho
mds1_cor_long %>%
  filter(Var1 == "x" & Var2 != "x") %>%
  arrange(-abs(cor)) %>%
  head(n_features) %>%
  mutate(cor = round(cor, 3))

features_y <- mds1_cor_long %>%
  filter(Var1 == "y" & Var2 != "y") %>%
  arrange(-abs(cor)) %>%
  head(n_features) %>%
  pull(Var2) %>%
  as.character()

# plot taxa vs mds1
global_F_long <- melt(data.matrix(global_F_rel_filt))
names(global_F_long) <- c("feature", "sample", "rel_abundance")
mds_taxa <- merge(mds_meta, global_F_long, by = "sample")

mds_taxa_plot <- mds_taxa %>%
  filter(feature %in% features_x)
mds_taxa_plot$feature <- factor(mds_taxa_plot$feature, levels = features_x)

scatter_F <- ggplot(mds_taxa_plot, aes(x = x, y = rel_abundance, color = site2)) +
  geom_point(size = 1) +
  facet_wrap(feature ~ ., scales = "free", ncol = 1, strip.position = "left") +
  scale_color_manual(values = global_pal) +
  theme_cowplot(11) +
  theme(axis.title = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, face = "italic", size = 9),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank()) +
  background_grid(major = "x")

# panel C: mds axes ----
mds1 <- ggplot(mds_meta, aes(site2, x, fill = site2)) +
  geom_jitter(alpha = 0.75, color = "darkgray", width = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  scale_fill_manual(values = global_pal) +
  theme_cowplot(12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  labs(x = "",
       y = "MDS 1") +
  background_grid(major = "y")

mds2 <- ggplot(mds_meta, aes(site2, y, fill = site2)) +
  geom_jitter(alpha = 0.75, color = "darkgray", width = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  scale_fill_manual(values = global_pal) +
  theme_cowplot(12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  labs(x = "",
       y = "MDS 2") +
  background_grid(major = "y")

# panel D: global shannon diversity ----

# shannon diversity
global_S_rare <- rrarefy(t(global_S), min(colSums(global_S)))

shannon_div <- diversity(global_S_rare, index = "shannon")
div <- data.frame("shannon_div" = shannon_div, "sample" = names(shannon_div))
div_meta <- merge(div, pheno_global, by = "sample")

div_meta <- div_meta[complete.cases(div_meta), ]

# add p-values
pvals <- compare_means(shannon_div ~ site2, data = div_meta, method = "wilcox.test", p.adjust.method = "fdr")

div_meta <- div_meta[complete.cases(div_meta), ]

shannon_global <- ggplot(div_meta, aes(site2, shannon_div, fill = site2)) + 
  geom_jitter(alpha = 0.75, color = "darkgray", width = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  labs(x = "",
       y = "Shannon Diversity",
       fill = "") +
  scale_fill_manual(values = global_pal) +
  scale_y_continuous(breaks = seq(1, 6, 1)) +
  theme_cowplot(12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  background_grid(major = "y")

# plot panels ----
col1 <- plot_grid(
  cohort_plot + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")),
  mds1 + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")),
  mds2 + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")),
  shannon_global + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")),
  ncol = 1,
  align = "v",
  axis = "lbr",
  labels = c("A", "C", "", "D"),
  rel_heights = c(0.3, 0.2, 0.2, 0.3)
)

mds_scatter <- plot_grid(
  global_mds + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"),
                     legend.position = "none"),
  scatter_F + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"),
                    legend.position = "none"),
  align = "hv",
  axis = "lr",
  ncol = 1,
  labels = c("B", ""),
  rel_heights = c(0.5, 0.5)
)

p1 <- plot_grid(
  col1,
  mds_scatter,
  nrow = 1,
  rel_widths = c(0.45, 0.55)
)

plot_grid(get_legend(global_mds), p1, ncol = 1, rel_heights = c(0.05, 0.95))

ggsave(here("final_plots/figure_3.png"), width = 9, height = 10, bg = "white")

