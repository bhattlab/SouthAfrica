library(gtools)
library(here)
library(tidyverse)

##### EDIT HERE: #####
ega_path <- ""
######################

ega_data <- read.csv(ega_path)

study_samples <- read.table(here("input_final/pheno/study_samples.txt"),
                            sep = "\t", header = T)

za_labels <- read.table(here("input_final/pheno/za_labels.tsv"), sep = "\t",
                        header = T)

za_labels[which(za_labels$sample == "A15"), "study_id"] <- "2554637a"

za_pheno <- study_samples %>%
  filter(keep) %>%
  left_join(ega_data, by = c("study_id" = "SampleID")) %>%
  select(sample, study_id, site, mb_age:glucose_reading) %>%
  left_join(za_labels, by = c("sample", "study_id", "site"))

za_meta <- za_pheno %>%
  select(sample, study_id, site)

save(za_meta, za_pheno, file = here("RData/metadata.RData"))
