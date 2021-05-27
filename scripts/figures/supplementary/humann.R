library(cowplot)
library(genefilter)
library(ggpubr)
library(here)
library(tidyverse)

## load data ----
source(here("scripts/load_data.R"))

## humann output ----

humann <- read.table(here("input_final/humann/join_metacyc-pwy-CPM.tsv"),
                     sep = "\t", comment.char = "", header = T)

names(humann) <- gsub("_concat_Abundance.CPM", "", names(humann))

humann <- humann %>%
  filter(!grepl("UNMAPPED|UNINTEGRATED|\\|", X..Pathway)) %>%
  column_to_rownames("X..Pathway")

humann <- humann[, za_pheno$sample]

# prevalence filter
humann_filt <- humann[genefilter(humann, pOverA(0.2, 0)), ]

humann_long <- humann_filt %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(cols = -Pathway, names_to = "sample") %>%
  merge(za_pheno, by = "sample")

wilcox_res <- humann_long %>%
  group_by(Pathway) %>%
  summarize(pvalue = wilcox.test(value ~ site)$p.value,
            pvalue2 = wilcox.test(value ~ site, alternative = "l")$p.value) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pvalue, method = "fdr")) %>%
  arrange(pvalue) %>%
  mutate(across(where(is.numeric), ~signif(., 2)))

