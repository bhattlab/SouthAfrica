# Analyses and plots for South African microbiome pilot data
# Fiona Tamburini

library(cowplot)
library(data.table)
library(genefilter)
library(tidyverse)

# setwd("/Users/tamburif/za-microbiome/")
setwd("/labs/asbhatt/fiona/za/za-microbiome/")

# metadata ----
za_meta <- readRDS("rds/za_meta.rds")
za_pheno <- readRDS("rds/za_pheno.rds")

# functions ----

# wilcoxon test function
humann_wilcox <- function(db, humann_path = "input_final/humann/full/join_uniref90_"){
  
  # test
  # db <- "rxn"
  
  humann_out <- data.frame(fread(paste0(humann_path, db, ".tsv"), sep = "\t", quote = "", header = T))
  
  names(humann_out) <- gsub("_concat_Abundance.*", "", names(humann_out))
  
  humann_out <- humann_out %>%
    filter(!grepl("UNMAPPED|UNGROUPED|\\|", X..Gene.Family)) %>%
    column_to_rownames("X..Gene.Family")
  
  humann_out <- humann_out[, za_meta$sample]
  
  # filter on abundance, prevalence
  humann_filt <- humann_out[genefilter(humann_out, pOverA(0.2, 0)), ]
  
  # variance filter
  humann_filt <- data.frame(varFilter(as.matrix(humann_filt), var.cutoff = 0.25))
  
  humann_filt <- humann_filt[, which(colSums(humann_filt) > 0)]
  
  h_long <- reshape2::melt(data.matrix(humann_filt))
  
  names(h_long) <- c("fxn", "sample", "value")
  
  h_long <- merge(h_long, za_meta, all.x = T, all.y = F, by = "sample")  
  
  wilcox_res <- h_long %>%
  group_by(fxn) %>%
  summarize(pvalue = wilcox.test(value ~ site)$p.value,
            mean_swt = mean(value[site == "Soweto"]),
            mean_bbr = mean(value[site == "Bushbuckridge"])) %>%
  mutate(padj = p.adjust(pvalue, method = "bonferroni"))

  # res_signif <- res %>% filter(p.adj < 0.05)
}

# run on all humann outputs ----
for (db in c("go", "ko", "eggnog", "pfam", "rxn", "level4ec")){
  
  res <- humann_wilcox(db)
  
  saveRDS(res, paste0("rds/humann/", db, "_res.rds"))
}

