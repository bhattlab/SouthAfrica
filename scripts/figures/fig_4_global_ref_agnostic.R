library(cowplot)
library(ggpubr)
library(harrietr)
library(here)
library(MASS)
library(reshape2)
library(tidyverse)
library(vegan)

## load data ----
source(here("scripts/load_data.R"))

## A) classification plot ----
# read data with unclassified reads
unclass <- read.table(here("input_final/readcounts/global_unclassified.txt"),
                      sep = "\t")
names(unclass) <- c("sample", "unclassified_reads")

# rampelli
unclass_rampelli <- read.table(
  here("input_final/readcounts/global_unclassified_rampelli.txt"),
  sep = "\t")
names(unclass_rampelli) <- c("sample", "unclassified_reads")

unclass <- rbind(unclass, unclass_rampelli)

readcounts <- read.table(here("input_final/readcounts/readcounts_global.txt"),
                         sep = "\t")
names(readcounts) <- c("sample", "total_reads")

readcounts_rampelli <- read.table(
  here("input_final/readcounts/readcounts_global_rampelli.txt"), sep = "\t")
names(readcounts_rampelli) <- c("sample", "total_reads")

readcounts <- rbind(readcounts, readcounts_rampelli)

# kracken/bracken results are in read pairs
unclass$unclassified_reads <- unclass$unclassified_reads * 2

unclass$total_reads <- readcounts[match(unclass$sample, readcounts$sample), "total_reads"]

unclass <- unclass %>%
  mutate(
    classified_reads = total_reads - unclassified_reads,
    classified_perc = classified_reads / total_reads
  )

# merge pheno data
unclass_pheno <- merge(pheno_global, unclass, by = "sample")

## Plot by site with statistical tests
# https://github.com/kassambara/ggpubr/wiki/Adding-Adjusted-P-values-to-a-GGPlot

my_comparisons <- list(
  c("Bushbuckridge", "Soweto"),
  # c("Bushbuckridge", "United States"),
  c("Soweto", "United States")
)

pvals <- unclass_pheno %>%
  rstatix::pairwise_wilcox_test(classified_perc ~ site2, comparisons = my_comparisons) %>%
  rstatix::add_y_position()
pvals$y.position <- pvals$y.position + 2

read_class_plot <- ggplot(unclass_pheno, aes(site2, classified_perc * 100)) + 
  # geom_boxplot(outlier.shape = NA, aes(fill = site2)) +
  # geom_jitter(position = position_jitterdodge(jitter.width = 2),
  #             alpha = 0.5, pch = 21, aes(fill = site2)) +
  geom_jitter(alpha = 0.75, color = "darkgray", width = 0.3) +
  geom_boxplot(outlier.shape = NA, aes(fill = site2), alpha = 0.75) +
  labs(x = "",
       y = "Classified reads (%)",
       fill="Site") +
  scale_fill_manual(values = global_pal) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  theme_cowplot() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.margin = unit(c(0, 0.5, 0, 0.5), unit = "cm")
  ) +
  stat_pvalue_manual(pvals, label = "p.adj.signif", tip.length = 0.01,
                     vjust = 0.5, hide.ns = T, y.position = c(100, 104)) +
  background_grid(major = "y")


## B) k-mer comparisons ----
# perform kmer nmds/mds
kmer_mds <- function(k){
  k_matrix_filt <- get(paste0("sourmash_k", k, "_track_abund"))
  
  k_dist <- as.dist(1 - k_matrix_filt)
  
  cmdscale(k_dist, eig = TRUE, x.ret = TRUE)
  
  isoMDS(k_dist)
}

# plot kmer mds
plot_kmer_mds <- function(mds, k){
  
  mds <- data.frame(mds$points)
  mds$sample <- row.names(mds)
  mds_pheno <- merge(mds, pheno_global, by = "sample", all.x = T, all.y = F)
  
  # randomize sample order for more even plotting
  set.seed(42)
  rows <- sample(nrow(mds_pheno))
  mds_pheno <- mds_pheno[rows, ]
  
  kmer_mds <- ggplot(mds_pheno, aes(X1, X2, color = site2)) +
    geom_point(size = 1.5, alpha = 0.75) +
    scale_color_manual(values = global_pal, guide = guide_legend(nrow = 1)) +
    labs(
      # title = paste("k=", k, sep = ""),
      x = "NMDS 1",
      y = "NMDS 2",
      color="",
      fill="Site"
    ) +
    theme_cowplot() +
    theme(
      legend.position = "top",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), unit = "cm")
    ) +
    annotate("text", label = paste0("italic(k) ==", k),
             x = min(mds_pheno$X1) + 0.1,
             y = max(mds_pheno$X2), parse = T) +
    background_grid() +
    coord_fixed()
  
  kmer_mds +
    stat_ellipse(aes(color = site2), type = 't', size = 0.5, alpha = 0.75,
                 show.legend = F)
}

k21_mds <- kmer_mds("21")
k31_mds <- kmer_mds("31")
k51_mds <- kmer_mds("51")

k21 <- plot_kmer_mds(k21_mds, 21)
k31 <- plot_kmer_mds(k31_mds, 31)
k51 <- plot_kmer_mds(k51_mds, 51)


## C) compare k-mer and bray pairwise distances ----
# species bray-curtis vs kmer angular similarity
for (method in c("bray", "jaccard")){
  
  suffix <- ifelse(method == "bray", "_track_abund", "")
  
  # veg_dist <- vegdist(t(global_S_css), method = method)
  veg_dist <- vegdist(t(global_S_rel), method = method)
  dist_long <- melt_dist(as.matrix(veg_dist))
  dist_long <- merge(dist_long, pheno_global[, -2], by.x = "iso1",
                     by.y = "sample", all.x = T, all.y = F)
  dist_long <- merge(dist_long, pheno_global[, -2], by.x = "iso2",
                     by.y = "sample", all.x = T, all.y = F)
  pair_dis_sp <- filter(dist_long, site2.x == site2.y)
  pair_dis_sp$method <- "Species"
  
  k_long <- melt_dist(1 - get(paste0("sourmash_k31", suffix)))
  k_long <- merge(k_long, pheno_global[, -2], by.x = "iso1",
                  by.y = "sample", all.x = T, all.y = F)
  k_long <- merge(k_long, pheno_global[, -2], by.x = "iso2",
                  by.y = "sample", all.x = T, all.y = F)
  pair_dis_k <- filter(k_long, site2.x == site2.y)
  pair_dis_k$method <- "K-mer"
  
  pair_dis <- rbind(pair_dis_sp, pair_dis_k)
  pair_dis$method <- factor(pair_dis$method, levels = c("Species", "K-mer"))
  
  assign(paste0("pair_dis", suffix), pair_dis)
}

## plot all species / k-mer comparisons
my_comparisons <- list(
  c("Soweto", "Tanzania"),
  c("Soweto", "Madagascar"),
  c("Soweto", "Bushbuckridge"),
  c("Soweto", "Sweden"),
  c("Soweto", "United States")
)

dist_plot_all <- ggplot(pair_dis_track_abund,
                        aes(x = site2.y,y = dist, fill = site2.y)) +
  geom_boxplot(outlier.size = 0.75) +
  facet_wrap(~method, scales = "free",
             labeller = labeller(groupwrap = label_wrap_gen(10))) +
  scale_fill_manual(values = global_pal) +
  labs(
    x = "",
    y = "Pairwise Distance"
  ) +
  theme_cowplot(12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",
                     method = "wilcox.test", vjust = 0.75) +
  background_grid(major = "y")


## only za/Sweden
my_comparisons <- list(
  c("Bushbuckridge", "Sweden"),
  c("Soweto", "Sweden")
)

dat_sp <- pair_dis_track_abund %>%
  filter(site2.y %in% c("Bushbuckridge", "Sweden"), method == "Species")
wilcox.test(dist ~ site2.y, dat_sp, alternative = "l")

dat_k <- pair_dis_track_abund %>%
  filter(site2.y %in% c("Bushbuckridge", "Sweden"), method == "K-mer")
wilcox.test(dist ~ site2.y, dat_k, alternative = "g")

## are dist distributions normal?
# d <- dat %>% filter(site2.y == "Bushbuckridge") %>% pull(dist)
# ggdensity(d)
# ggqqplot(d)
# shapiro.test(sample(d, 5000))

za_sw <- pair_dis_track_abund %>%
  filter(site2.y %in% c("Bushbuckridge", "Sweden"))

dist_plot_za_swed <- ggplot(za_sw, aes(x = site2.y, y = dist, fill = site2.y)) +
  geom_boxplot(outlier.size = 0.75) +
  facet_wrap(~method, scales = "free",
             labeller = labeller(groupwrap = label_wrap_gen(10))) +
  scale_fill_manual(values = global_pal[c(3, 5)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.075))) +
  # expand_limits(y = 1.1) +
  labs(
    x = "",
    y = "Pairwise Distance"
  ) +
  theme_cowplot(12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  # stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test") +
  stat_compare_means(comparisons = list(c("Bushbuckridge", "Sweden")),
                     label = "p.signif", method = "wilcox.test", vjust = 0.5) +
  background_grid(major = "y")


## Madagascar/USA
mad_usa <- pair_dis_track_abund %>%
  filter(site2.y %in% c("Madagascar", "United States"))

dist_plot_mad_usa <- ggplot(mad_usa, aes(x = site2.y, y = dist, fill = site2.y)) +
  geom_boxplot(outlier.size = 0.75) +
  facet_wrap(~ method, scales = "free",
             labeller = labeller(method = label_wrap_gen(10))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.075))) +
  scale_fill_manual(values = global_pal[c(2, 6)]) +
  labs(
    x = "",
    y = "Pairwise Distance"
  ) +
  theme_cowplot(12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  stat_compare_means(comparisons = list(c("Madagascar", "United States")),
                     label = "p.signif", method = "wilcox.test", vjust = 0.5) +
  # geom_signif(comparisons = list(c("Madagascar", "United States")),
  #             map_signif_level = TRUE, test = "wilcox.test") +
  background_grid(major = "y")

row1 <- plot_grid(get_legend(k21))

row2 <- plot_grid(
  read_class_plot,
  k31 + theme(legend.position = "none"),
  labels = c("A", "B"),
  nrow = 1
)

row4 <- plot_grid(
  dist_plot_za_swed,
  dist_plot_mad_usa,
  labels = c("D", "E"),
  align = "hv",
  axis = "bt"
)

plot_grid(
  row1,
  row2,
  dist_plot_all,
  row4,
  labels = c("", "", "C", ""),
  ncol = 1,
  rel_heights = c(0.05, 0.3, 0.35, 0.3)
)

ggsave(here("final_plots/figure_4.png"), width = 8, height = 12.5)


### supplementary: Jaccard ----

## plot all species / k-mer comparisons
ggplot(pair_dis, aes(x = site2.y, y = dist, fill = site2.y)) +
  geom_boxplot(outlier.size = 0.75) +
  facet_wrap(~method, scales = "free", labeller = labeller(groupwrap = label_wrap_gen(10))) +
  scale_fill_manual(values = global_pal) +
  labs(
    x = "",
    y = "Pairwise Jaccard Distance"
  ) +
  theme_cowplot(12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  background_grid(major = "y")

ggsave(here("final_plots/misc/jaccard_dist_supp.png"), width = 8, height = 5)
