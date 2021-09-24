library(genefilter)
library(gtools)
library(here)
library(tidyverse)

## supplementary feature tables ----

## load data ----
load(here("RData/metadata.RData"))
load(here("RData/sequence_data.RData"))

# only include features present at at least 0.01% relab in at least one sample
labs <- za_pheno %>%
  select(sample, label) %>%
  as.data.frame

# genus
genus_supp <- za_G_rel * 100
names(genus_supp) <- labs[match(names(genus_supp), labs$sample), "label"]
genus_supp <- genus_supp[, mixedorder(names(genus_supp))]

genus_supp <- genus_supp[genefilter(genus_supp, kOverA(1, 0.01)), ]

genus_supp <- genus_supp %>%
  rownames_to_column("Genus")

write.table(genus_supp, "final_tables/table_s4_genus_classification_rel.txt",
            sep = "\t", row.names = F, quote = F)

# species
species_supp <- za_S_rel * 100
names(species_supp) <- labs[match(names(species_supp), labs$sample), "label"]
species_supp <- species_supp[, mixedorder(names(species_supp))]

species_supp <- species_supp[genefilter(species_supp, kOverA(1, 0.01)), ]

species_supp <- species_supp %>%
  rownames_to_column("Species")

write.table(species_supp, "final_tables/table_s5_species_classification_rel.txt",
            sep = "\t", row.names = F, quote = F)
