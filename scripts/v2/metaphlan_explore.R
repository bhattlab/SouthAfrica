library(tidyverse)
library(harrietr)

setwd("/Users/tamburif/za-microbiome")

# read metadata ---

za_meta <- readRDS("rds/za_meta.rds")
za_pheno <- readRDS("rds/za_pheno.rds")

# read metaphlan data ---
m_out <- read.table("input_final/metaphlan_merge.txt", sep = "\t", header = T)

colnames(m_out) <- gsub("_concat.+", "", colnames(m_out))

m_out_s <- m_out %>%
  filter(grepl("s__", clade_name)) %>%
  dplyr::select(-NCBI_tax_id) %>%
  column_to_rownames("clade_name")

# keep relevant samples
m_out_s <- m_out_s[, za_meta$sample]

# rm zero features
m_out_s <- m_out_s[rowSums(m_out_s) > 0, ]

# most abundant species across population ----
n <- 10
sp <- names(rev(sort(rowMeans(m_out_s)))[1:n])
# write.table(sp, "/Users/tamburif/Desktop/sp.txt", sep = "\t", quote = F, row.names = F)

# wilcoxon test for diff bugs across site ----

m_long <- reshape2::melt(data.matrix(m_out_s))

names(m_long) <- c("species", "sample", "rel_abund")

m_long <- merge(m_long, za_pheno, all.x = T, all.y = F, by = "sample")  

# wilcox_res <- m_long %>%
#   group_by(species) %>%
#   do(w = wilcox.test(rel_abund ~ site, data = ., paired = FALSE)) %>%
#   summarise(species, pvalue = w$p.value) %>%
#   mutate(padj = p.adjust(pvalue, method = "fdr"))

wilcox_res <- m_long %>%
  group_by(species) %>%
  summarize(pvalue = wilcox.test(rel_abund ~ site)$p.value,
            mean_swt = mean(rel_abund[site == "Soweto"]),
            mean_bbr = mean(rel_abund[site == "Bushbuckridge"])) %>%
  mutate(padj = p.adjust(pvalue, method = "fdr"))

wilcox_res %>%
  arrange(padj) %>%
  filter(padj < 0.1) %>%
  mutate(sname = gsub(".+\\|", "", species),
         inc_bbr = mean_bbr > mean_swt)

# continuous associations with bmi ----

# m_out_pseudo <- m_out_s
# m_out_pseudo[m_out_pseudo == 0] <- 1e-4
# m_t <- data.frame(t(log2(m_out_s)))

m_t <- data.frame(t(m_out_s[, !names(m_out_s) == "C9"]))
m_t$bmi <- za_pheno[match(rownames(m_t), za_pheno$sample), "bmi"]

scor <- cor(m_t, method = "spearman")

scor_long <- melt_dist(scor)

scor_long %>%
  filter(iso1 == "bmi" | iso2 == "bmi", abs(dist) > 0.2) %>%
  arrange(-dist) %>%
  mutate(iso1 = gsub(".+\\.", "", iso1),
         iso2 = gsub(".+\\.", "", iso2))

ggplot(m_t, aes(bmi, k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfovibrionales.f__Desulfovibrionaceae.g__Desulfovibrio.s__Desulfovibrio_piger)) +
  geom_point()

# categorical associations with lean/obese ----




