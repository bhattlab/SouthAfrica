library(cowplot)
library(genefilter)
library(ggpubr)
library(harrietr)
library(here)
library(tidyverse)
library(Maaslin2)

# load data ----
source(here("scripts/load_data.R"))

# read card metadata ----
aro_desc <- read.csv(here("input_final/shortbred/aro.csv"))

# read card shortbred data ---
card <- read.table(here("input_final/shortbred/CARD_merge.txt"),
                   sep = "\t", header = T)

card_wide <- card %>%
  filter(Sample %in% za_meta$sample, Count > 0, Hits > 0) %>%
  separate(col = Family, into = c("db", "accession", "aro", "name"),
           sep = "\\|") %>%
  select(aro, Sample, Count) %>%
  group_by(aro, Sample) %>%
  summarise(count = sum(Count)) %>%
  ungroup() %>%
  mutate(aro = gsub("_", ":", aro)) %>%
  pivot_wider(id_cols = aro, names_from = Sample, values_from = count,
              values_fill = 0) %>%
  left_join(aro_desc, by = c("aro" = "Accession")) %>%
  # shorten labels for plotting
  mutate(Name_full = Name,
         Name = ifelse(nchar(Name) > 30,
                       gsub("^(\\S+ \\S+ \\S+).+", "\\1\\.\\.\\.", Name),
                       Name),
         aro_desc = paste0(Name, " (", aro, ")"))

# table of shortened names
aro_names <- card_wide %>%
  filter(Name != Name_full) %>%
  select(aro, aro_desc, Name_full, Name)

card_wide <- card_wide %>%
  column_to_rownames("aro_desc") %>%
  select(one_of(za_meta$sample))

# remove low prevalence features
card_wide <- card_wide[genefilter(card_wide, pOverA(0.1, 0)), ]

aro_names <- aro_names %>%
  filter(aro_desc %in% rownames(card_wide)) %>%
  select(aro, Name_full)

# maaslin ----

# maaslin parameters
model <- "CPLM"
norm <- "TSS"
transf <- "NONE"

prefix <- paste("abx", model, norm, sep = "_")
dirname <- here("output/maaslin2", prefix)

meta <- za_meta
rownames(meta) <- meta$sample

card_maaslin <- card_wide[, meta$sample]

names_lookup <- card_maaslin %>%
  rownames_to_column("aro_desc") %>%
  select(aro_desc) %>%
  mutate(feature = make.names(aro_desc))

# run maaslin
fit_data <- Maaslin2(card_maaslin, meta,
                     dirname,
                     transform = transf,
                     analysis_method = model,
                     fixed_effects = c("site"),
                     normalization = norm,
                     standardize = FALSE,
                     plot_heatmap = F,
                     plot_scatter = F)

res <- fit_data$results %>%
  left_join(names_lookup, by = "feature") %>%
  as_tibble() %>%
  filter(metadata == "site") %>%
  arrange(pval) %>%
  mutate(model = model,
         norm = norm)

# maaslin coef plot ----

res_plot <- res %>%
  filter(qval < 0.05,
         abs(coef) > 0.5) %>%
  mutate(effect = factor(coef < 0),
         aro_desc = fct_reorder(aro_desc, coef))

a <- ggplot(res_plot, aes(coef, aro_desc, fill = effect)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rev(as.vector(za_pal)),
                    labels = c("Enriched in SWT", "Enriched in BBR")) +
  theme_cowplot() +
  background_grid() +
  labs(x = "Coefficient",
       y = "",
       fill = "") +
  theme(legend.position = "top",
        legend.justification = "center")

# card heatmap ----

set.seed(1)

card_mat <- card_wide

names(card_mat) <- za_pheno[match(names(card_mat), za_pheno$sample), ] %>%
  pull(label_abbrev)

# cluster samples
ddist <- dist(as.matrix(t(card_mat)), method = "canberra")
hc_samples <- hclust(d = ddist)

# cluster features
fdist <- dist(as.matrix(card_mat), method = "euclidean")
hc_features <- hclust(d = fdist)

card_long_pc <- card_wide %>%
  rownames_to_column("aro_desc") %>%
  pivot_longer(cols = -aro_desc, names_to = "sample", values_to = "rpkm") %>%
  left_join(za_pheno, by = "sample") %>%
  filter(aro_desc %in% rownames(card_mat)) %>%
  mutate(rpkm = ifelse(rpkm == 0, min(rpkm[rpkm > 0])/2, rpkm),
         label_abbrev = factor(label_abbrev,
                               levels = names(card_mat)[hc_samples$order]),
         aro_desc = factor(aro_desc,
                         levels = rownames(card_mat)[hc_features$order]))

b2 <- ggplot(card_long_pc, aes(label_abbrev, aro_desc, fill = log2(rpkm))) +
  geom_tile() +
  viridis::scale_fill_viridis() +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 6),
        plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm"),
        legend.position = "bottom",
        legend.justification = "center") +
  labs(x = "Participant",
       y = "")

b1 <- ggplot(card_long_pc, aes(label_abbrev, "1", fill = site)) +
  geom_tile() +
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"))

b <- plot_grid(b1, b2, ncol = 1, align = "v", axis = "lr",
               rel_heights = c(0.05, 0.95))

# table of labels that were truncated for plotting
lvls <- levels(card_long_pc$aro_desc)
lvls <- lvls[grepl("\\.\\.\\.", lvls)]
lvls <- rev(gsub(".+\\.\\.\\. \\((.+)\\)", "\\1", lvls))

aro_tbl <- aro_names %>%
  mutate(aro = factor(aro, levels = lvls)) %>%
  arrange(aro) %>%
  ggtexttable(rows = NULL,
              theme = ttheme("classic", base_size = 10,
                             padding = unit(c(4, 4), "mm")),
              cols = c("Accession", "Name"))

# plot_grid(b, a, ncol = 1, rel_heights = c(0.6, 0.4), labels = c("A", "B"))
p2 <- plot_grid(a, aro_tbl, rel_widths = c(0.55, 0.45), labels = c("B", "C"))

plot_grid(b, p2, ncol = 1, rel_heights = c(0.6, 0.4), labels = c("A", ""))

ggsave(here("final_plots/supplementary/figure_S8_card.png"),
       width = 16, height = 18, bg = "white")
