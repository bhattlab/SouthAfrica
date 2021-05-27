## update Voke metadata table for EGA submission

library(dplyr)

setwd("/Users/tamburif/za-microbiome/misc/ENA metadata")

ega_meta <- read.csv("agt_swt_metadata_study_rev_1.csv")

# add patient JXTVM
za_pheno <- read.table("/Users/tamburif/za-microbiome/input_final/za_pheno.txt", sep = "\t", header = T, quote = "")

meta_add <- za_pheno %>%
  filter(study_id == "JXTVM") %>%
  mutate(
    batch = "b2",
    bmi = NA,
    bmi_info = NA
    ) %>%
  select(study_id, batch, bmi, bmi_info, hiv_result, glucose_reading, site)

names(meta_add) <- names(ega_meta)

ega_meta <- rbind(ega_meta, meta_add)


# add blood pressure measurements
za_pheno_full <- read.table("/Users/tamburif/za-microbiome/misc/ENA metadata/za_pheno_all.txt", sep = "\t", header = T, quote = "")

# fix sample ID
ega_meta[ega_meta$SampleID == "2554637a", "SampleID"] <- "2554637"

ega_meta$SampleID[!(ega_meta$SampleID %in% za_pheno_full$study_id)]

ega_meta_bp <- merge(ega_meta, za_pheno_full %>% select(study_id, starts_with("bp")), by.x = "SampleID", by.y = "study_id", all.x = T, all.y = F)

ega_meta_bp[ega_meta_bp$SampleID == "2554637", "SampleID"] <- "2554637a"

write.csv(ega_meta_bp, "/Users/tamburif/za-microbiome/misc/ENA metadata/agt_swt_metadata_study_rev_1_FT.csv", quote = F)