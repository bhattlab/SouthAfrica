library(cowplot)
library(genefilter)
library(ggpubr)
library(harrietr)
library(here)
library(tidyverse)

## load data ----
source(here("scripts/load_data.R"))

## read vfdb metadata ----
vfdb_meta <- read.table(here("input_final/shortbred/vfdb_meta.txt"),
                        sep = "\t", quote = "", col.names = c("Family", "gene", "description"))
victors_meta <- read.table(here("input_final/shortbred/victors_meta.txt"),
                           sep = "\t", comment.char = "", quote = "")
victors_meta <- victors_meta %>%
  mutate(Family = paste(V1, V2, V3, V4, sep = "|")) %>%
  dplyr::select(Family, V5)
names(victors_meta) <- c("Family", "description")

## read shortbred data ---
card <- read.table(here("input_final/shortbred/CARD_merge.txt"), sep = "\t", header = T)
vf <- read.table(here("input_final/shortbred/VF_merge_filtered.txt"), sep = "\t", header = T)

card_wide <- card %>%
  filter(Sample %in% za_meta$sample, Count > 0, Hits > 0) %>%
  pivot_wider(id_cols = c(Family, TotMarkerLength), names_from = Sample, values_from = Count, values_fill = 0) %>%
  separate(col = Family, into = c("db", "accession", "aro", "name"), sep = "\\|", remove = F)

vf_wide <- vf %>%
  filter(Sample %in% za_meta$sample) %>%
  pivot_wider(id_cols = c(Family, TotMarkerLength), names_from = Sample, values_from = Count, values_fill = 0)

## filter low prevalence
card_wide <- card_wide[genefilter(card_wide[, 7:ncol(card_wide)], pOverA(0.1, 0)), ]
vf_wide <- vf_wide[genefilter(vf_wide[, 3:ncol(vf_wide)], pOverA(0.1, 0)), ]

# remove mvirdb - no metadata
vf_wide <- vf_wide %>%
  filter(!grepl("^virulence\\|", Family))

# VF labels
meta_both <- rbind(vfdb_meta[, -2], victors_meta)
meta_both$Family <- gsub(":", "_", meta_both$Family)

vf_wide <- vf_wide %>%
  arrange(Family) %>%
  mutate(Family_label = row_number()) %>%
  relocate(Family, Family_label, TotMarkerLength)

vf_tbl <- vf_wide %>%
  select(Family_label, Family) %>%
  distinct() %>%
  mutate(db = gsub("\\|.+", "", Family),
         Family_match = gsub("VFDB\\||Victors\\||\\|$", "", Family), # format for matching to meta
         Family_match = gsub("_([0-9]$)", "\\.\\1", Family_match), # format to match meta
         Family_match = gsub("(_virulence|_Victors|_gi).+", " ", Family_match), # format to handle multiple gi ids in one entry
         Family_match = gsub("|\\s+$|\\|$", "", Family_match),
         description = meta_both[match(Family_match, meta_both$Family), "description"]) %>%
  arrange(Family_label) %>%
  select(Family_label, Family, description)
names(vf_tbl) <- c("VF", "Shortbred VF family", "Description")

write.table(vf_tbl, here("final_tables/shortbred_VF.txt"), sep = "\t",
            quote = F, row.names = F)

## read card metadata ----
aro_desc <- read.csv(here("input_final/shortbred/aro.csv"))
aro_desc$Accession <- gsub(":", "_", aro_desc$Accession)

## card wilcoxon tests ----

# by site
card_long <- reshape2::melt(card_wide,
                            id.vars = c("Family", "db", "accession", "aro", "name", "TotMarkerLength"),
                            variable.name = "sample", value.name = "rpkm")

card_long <- merge(card_long, za_pheno, all.x = T, all.y = F, by = "sample")  

wilcox_res <- card_long %>%
  group_by(Family, aro, name) %>%
  summarize(pvalue = wilcox.test(rpkm ~ site)$p.value,
            mean_swt = mean(rpkm[site == "Soweto"]),
            mean_bbr = mean(rpkm[site == "Bushbuckridge"])) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pvalue, method = "fdr"),
         aro_name = aro_desc[match(aro, aro_desc$Accession), "Name"],
         inc_bbr = mean_bbr > mean_swt)

# obese vs lean
wilcox_res_bmi <- card_long %>%
  filter(bmi_status != "Overweight") %>%
  group_by(Family, aro, name) %>%
  summarize(pvalue = wilcox.test(rpkm ~ bmi_status)$p.value,
            mean_obese = mean(rpkm[bmi_status == "Obese"]),
            mean_lean = mean(rpkm[bmi_status == "Lean"])) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pvalue, method = "fdr"),
         aro_name = aro_desc[match(aro, aro_desc$Accession), "Name"],
         inc_obese = mean_obese > mean_lean)

## plot card ----
fdr_cutoff <- 0.01
res_top <- wilcox_res %>%
  filter(padj < fdr_cutoff) %>%
  arrange(inc_bbr, padj)

card_plot <- card_long %>%
  filter(Family %in% res_top$Family) %>%
  mutate(label = paste(aro, name, sep = "\n"))

# pseudo count
card_plot$rpkm[card_plot$rpkm == 0] <- 0.01

a <- ggplot(card_plot, aes(site_abbrev, rpkm, fill = site)) +
  geom_jitter(width = 0.3, alpha = 0.85, size = 1, color = "gray50") +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  facet_wrap(~label, scales = "free_y", ncol = 8) +
  scale_y_log10(label = scales::comma) +
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 11),
        legend.position = "none") +
  background_grid(major = "y") +
  labs(x = "Site",
       y = "RPKM")

wilcox_res %>%
  mutate(mean_bbr_p = ifelse(mean_bbr == 0, 0.0001, mean_bbr),
         mean_swt_p = ifelse(mean_swt == 0, 0.0001, mean_swt),
         log2fc = log2(mean_swt/mean_bbr)) %>%
  arrange(-abs(log2fc)) %>%
  filter(padj < 0.05)

## heatmap
card_mat <- card_wide %>%
  column_to_rownames("Family") %>%
  select(one_of(za_meta$sample))

names(card_mat) <- za_pheno[match(names(card_mat), za_pheno$sample), ] %>%
  pull(label_abbrev)

# keep features >0 rpkm in 10% of individuals
# card_mat <- card_mat[genefilter(card_mat, pOverA(0.1, 0)), ]

# cluster samples
ddist <- dist(as.matrix(t(card_mat)), method = "canberra")
hc_samples <- hclust(d = ddist)

# cluster features
fdist <- dist(as.matrix(card_mat), method = "euclidean")
hc_features <- hclust(d = fdist)

card_long_pseudo <- card_long %>%
  filter(Family %in% rownames(card_mat))

card_long_pseudo$rpkm[card_long_pseudo$rpkm == 0] <- 0.01

card_long_pseudo$label_abbrev <- factor(card_long_pseudo$label_abbrev,
                                  levels = names(card_mat)[hc_samples$order])

card_long_pseudo$Family <- factor(card_long_pseudo$Family,
                                  levels = rownames(card_mat)[hc_features$order])


b2 <- ggplot(card_long_pseudo, aes(label_abbrev, Family, fill = log2(rpkm))) +
  geom_tile() +
  viridis::scale_fill_viridis() +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
  labs(x = "Participant",
       y = "")

b1 <- ggplot(card_long_pseudo, aes(label_abbrev, "1", fill = site)) +
  geom_tile() +
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = "none")

b <- plot_grid(b1, b2, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.05, 0.95))

plot_grid(b, a, ncol = 1, rel_heights = c(0.6, 0.4), labels = c("A", "B"))

ggsave(here("final_plots/supplementary/figure_S8_shortbred_card.png"),
       width = 16, height = 18)

# clst <- c("C14", "B11", "B16", "B33", "C1", "C19", "A69", "B9", "B20", "B26")

## vfdb ----
### VF by site ----
vf_long <- reshape2::melt(vf_wide,
                          id.vars = c("Family", "Family_label", "TotMarkerLength"),
                          variable.name = "sample", value.name = "rpkm")

vf_long <- merge(vf_long, za_pheno, all.x = T, all.y = F, by = "sample")  

wilcox_res <- vf_long %>%
  group_by(Family) %>%
  summarize(pvalue = wilcox.test(rpkm ~ site)$p.value,
            mean_swt = mean(rpkm[site == "Soweto"]),
            mean_bbr = mean(rpkm[site == "Bushbuckridge"])) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pvalue, method = "fdr"))

fdr_cutoff <- 0.01
res_top <- wilcox_res %>%
  filter(padj < fdr_cutoff) %>%
  mutate(inc_bbr = mean_bbr > mean_swt) %>%
  arrange(inc_bbr, padj)

n <- 20
vf_plot <- vf_long %>%
  filter(Family %in% res_top$Family[1:n]) %>%
  mutate(label = as.numeric(as.factor(Family)))

# pseudo count
vf_plot$rpkm[vf_plot$rpkm == 0] <- 0.01

a <- ggplot(vf_plot, aes(site_abbrev, rpkm, fill = site)) +
  geom_jitter(width = 0.3, alpha = 0.85, size = 1, color = "gray50") +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  facet_wrap(~label, scales = "free_y", ncol = 9) +
  scale_y_log10(label = scales::comma) +
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 11),
        legend.position = "none") +
  background_grid(major = "y") +
  labs(x = "Site",
       y = "RPKM")

## heatmap
vf_mat <- vf_wide %>%
  column_to_rownames("Family") %>%
  select(one_of(za_meta$sample))

names(vf_mat) <- za_pheno[match(names(vf_mat), za_pheno$sample), ] %>%
  pull(label_abbrev)

# keep features >0 rpkm in 10% of individuals
vf_mat <- vf_mat[genefilter(vf_mat, pOverA(0.1, 0)), ]

# cluster samples
ddist <- dist(as.matrix(t(vf_mat)), method = "canberra")
hc_samples <- hclust(d = ddist)

# cluster features
fdist <- dist(as.matrix(vf_mat), method = "euclidean")
hc_features <- hclust(d = fdist)

vf_long_pseudo <- vf_long %>%
  filter(Family %in% rownames(vf_mat))

vf_long_pseudo$rpkm[vf_long_pseudo$rpkm == 0] <- 0.01

# vf_long_pseudo$sample <- factor(vf_long_pseudo$sample,
#                                   levels = names(vf_mat)[hc_samples$order])
vf_long_pseudo$label_abbrev <- factor(vf_long_pseudo$label_abbrev,
                                levels = names(vf_mat)[hc_samples$order])

vf_long_pseudo$Family <- factor(vf_long_pseudo$Family,
                                  levels = rownames(vf_mat)[hc_features$order])


b2 <- ggplot(vf_long_pseudo, aes(label_abbrev, factor(Family_label), fill = log2(rpkm))) +
  geom_tile() +
  viridis::scale_fill_viridis() +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
  labs(x = "Participant",
       y = "")

b1 <- ggplot(vf_long_pseudo, aes(label_abbrev, "1", fill = site)) +
  geom_tile() +
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = "none")


b <- plot_grid(b1, b2, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.05, 0.95))

plot_grid(b, a, ncol = 1, rel_heights = c(0.8, 0.2), labels = c("A", "B"))

ggsave(here("final_plots/supplementary/shortbred_vf.png"),
       width = 16, height = 18)

