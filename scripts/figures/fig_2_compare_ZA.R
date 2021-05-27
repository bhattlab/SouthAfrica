library(cowplot)
library(DESeq2)
library(genefilter)
library(ggpubr)
library(here)
library(reshape2)
library(tidyverse)
library(vegan)

## load data ----
source(here("scripts/load_data.R"))

## A) MDS ordination ----
# find bray-curtis distance with vegdist
# vare_dis <- vegdist(t(za_S_rel), method = "bray")
vare_dis <- vegdist(t(za_S_css), method = "bray")

# adonis
meta <- za_meta %>% filter(sample %in% names(za_S_rel))
meta <- meta[match(names(za_S_rel), meta$sample), ]
adonis(vare_dis ~ site, data = meta)

dispersion <- betadisper(vare_dis, group = meta$site)
permutest(dispersion)

# calculate mds
mds <- cmdscale(vare_dis, eig = TRUE, x.ret = TRUE)

mds_values <- mds$points
wa_scores <- wascores(mds_values, t(za_S_rel))
wa_scores <- data.frame(sample = rownames(wa_scores),
                        x = wa_scores[,1],
                        y = wa_scores[,2])

# isolate taxa with strongest contribution to principal coordinate axes
n_taxa <- 10
wa_scores_1<- head(arrange(wa_scores, desc(abs(wa_scores$x))), n = n_taxa)
wa_scores_2<- head(arrange(wa_scores, desc(abs(wa_scores$y))), n = n_taxa)
wa_scores_final <- rbind(wa_scores_1, wa_scores_2)

# calculate percentage of variation that each MDS axis accounts for
mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)

# plot
mds_data <- data.frame(sample = rownames(mds_values),
                       x = mds_values[,1],
                       y = mds_values[,2])

# merge pheno data
mds_meta <- merge(mds_data, za_meta, by = "sample")

mds_plot <- ggplot(mds_meta, aes(x, y, color = site)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = za_pal) +
  labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
       y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
       color = "") +
  theme_cowplot() +
  background_grid() +
  coord_fixed() +
  theme(legend.position = "top",
        legend.justification = "center")

mds_plot_ci <- mds_plot +
  stat_ellipse(aes(color = site), type = 't', size = 1, show.legend = F)

## plot mds plot highlighting nanopore samples ----
nmag_samples <- c("C27", "C29", "C33")

mds_meta_nmag <- mds_meta %>%
  mutate(highlight = ifelse(sample %in% nmag_samples, T, F))

ggplot(mds_meta_nmag, aes(x, y, color = highlight, size = highlight, shape = site)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c("darkgrey", "red")) +
  scale_size_manual(values = c(2, 3)) +
  labs(
    x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
    y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
    shape = "Site"
  ) +
  theme_cowplot() +
  guides(size = F, color = F)

ggsave("final_plots/misc/nanopore_mds.png", width = 7, height = 5)


### B) Shannon diversity ----
# find shannon diversity with vegdist
shannon_div <- diversity(t(za_S_rel), index = "shannon")
div <- data.frame("shannon_div" = shannon_div, "sample" = names(shannon_div))
div_meta <- merge(div, za_meta, by = "sample")

# pval
# pval <- compare_means(shannon_div ~ site2, data = div_meta, method = "wilcox.test", p.adjust.method = "fdr")

div_plot <- ggplot(div_meta, aes(site, shannon_div)) + 
  # geom_boxplot(outlier.shape = NA, aes(fill = site)) +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.75), alpha = 0.5, aes(fill = site), pch = 21) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha = 0.85, aes(fill = site), color = "darkgray") +
  geom_boxplot(outlier.shape = NA, aes(fill = site), alpha = 0.5) +
  labs(x = "",
       y = "Shannon Diversity",
       fill="") +
  scale_color_manual(values = za_pal) +
  scale_fill_manual(values = za_pal) +
  scale_y_continuous(limits = c(2, 6.6)) +
  stat_compare_means(method = "wilcox", label = "p.signif",
                     comparisons = list(c("Soweto", "Bushbuckridge")), label.y = 6.5) +
  # geom_signif(comparisons = list(c("Soweto", "Bushbuckridge")), map_signif_level = TRUE, y_position = 6.5, test = "wilcox.test") +
  theme_cowplot(14) +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  ) +
  background_grid(major = "y")

wilcox.test(shannon_div ~ site, data = div_meta, exact = T)


## C) Differential features - DESeq2 ----
count_data <- za_G
count_data <- count_data[, za_meta$sample]

stopifnot(all(za_meta$sample == names(count_data)))

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = za_meta,
                              design = ~ site)

# estimate size factors
dds <- estimateSizeFactors(dds, type = "poscounts")

# filter normalized counts
idx <- genefilter(counts(dds, normalized = TRUE), pOverA(0.2, 500))
dds <- dds[idx, ]

# deseq2
dds <- DESeq(dds, parallel = T)

res <- results(dds, name = "site_Soweto_vs_Bushbuckridge", alpha = 0.05)

resOrdered <- res %>%
  as_tibble(rownames = "genus") %>%
  arrange(pvalue) %>%
  mutate(site = ifelse(log2FoldChange < 0, "Bushbuckridge", "Soweto"))

resOrdered %>%
  filter(abs(log2FoldChange) > 1, padj < 0.05) %>%
  arrange(log2FoldChange)

res_df_filt <- resOrdered %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(log2FoldChange)

res_df_filt$genus <- factor(res_df_filt$genus,
                            levels = as.character(res_df_filt$genus))

deseq_genera <- ggplot(res_df_filt, aes(log2FoldChange, genus, fill = site)) +
  # geom_point(size = 2, color = "black", pch = 21) +
  geom_bar(stat = "identity") +
  theme_cowplot() +
  theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = za_pal) +
  # scale_x_continuous(breaks = seq(-2, 2, 1)) +
  labs(
    x = "Log2 Fold Change (Soweto/Bushbuckridge)",
    y = "Genus",
    fill = ""
  ) +
  background_grid() +
  theme(legend.position = "top",
        legend.justification = "center")

res_table_G <- resOrdered %>%
  arrange(site, log2FoldChange) %>%
  mutate(across(where(is.numeric), ~ signif(., 3))) %>%
  filter(padj < 0.05) %>%
  select(genus, log2FoldChange, pvalue, padj) %>%
  arrange(-abs(log2FoldChange))

write.table(res_table_G, "final_tables/table_s7_deseq_genera.txt", sep = "\t",
            col.names = c("Genus", "Log2 fold change (SWT/BBR)", "P-value", "FDR-adjusted p-value"),
            row.names = F, quote = F)

## Figure 1
row1 <- plot_grid(mds_plot_ci, div_plot, labels = c('A', 'B'), align = "hv",
                  axis = "bt", label_size = 14, rel_widths = c(0.65, 0.35))

plot_grid(row1, deseq_genera, labels = c("", "C"), ncol = 1, label_size = 14,
          rel_heights = c(0.4, 0.6))

ggsave(here("final_plots/figure_2.png"), width = 9, height = 10)
