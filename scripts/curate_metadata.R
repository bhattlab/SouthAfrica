## Curate and combine microbiome pilot survey data from Agincourt batch 1 and 2, and Soweto collections

library(gtools)
library(here)
library(tidyverse)

# sequenced samples ----
sequenced_samples <- read.table(here("input_final/pheno/za_sequenced_samples.txt"),
                                sep = "\t", header = T)
sequenced_samples$sample <- gsub("(\\d+)([A-Z])", "\\2\\1", sequenced_samples$sample)

# remove duplicate
sequenced_samples <- sequenced_samples %>%
  filter(sample != "A14")

# voke metadata ----
voke_meta <- read.csv(here("input_final/pheno/agt_swt_metadata_study.csv"),
                      header = T)

# voke_meta %>%
#   group_by(site) %>%
#   tally()
# 
# voke_meta %>%
#   filter(is.na(hiv_result))

# from Voke email 9/10/18
agin_hiv <- c("MB17095", "MB17021", "MB17029", "MB17048", "MB17086", "MB17099",
              "MB17087", "MB17071", "MB17032", "MB17023")

# soweto pheno data ----
swt <- read.table(here("input_final/pheno/soweto.txt"), sep = "\t", header = T)

# swt %>%
#   group_by(HIV.test.results) %>%
#   tally()
# 
# swt %>%
#   filter(MicrobiomeStudy_ID %in% sequenced_samples$study_id) %>%
#   group_by(HIV.test.results) %>%
#   tally()

# bushbuckridge 1 ----
bbr_1 <- read.csv(here("input_final/pheno/agincourt-pilot-data.csv"), header = T)

# bbr_1 %>%
#   filter(!is.na(study_id) & study_id %in% sequenced_samples$study_id) %>%
#   group_by(hiv_pos) %>%
#   tally()
# 
# bbr_1 %>%
#   filter(!(study_id %in% sequenced_samples$study_id)) %>%
#   tally()

# bushbuckridge 2 ----
bbr_2 <- read.csv(here("input_final/pheno/redcap_batch2.csv"), header = T)

# bbr_2 %>%
#   filter(record_id %in% sequenced_samples$study_id) %>%
#   group_by(hiv_result) %>%
#   tally()

# merge metadata ----
# soweto
swt_merge <- swt %>%
  mutate(site = "Soweto") %>%
  select(-BTT_IDs)

names(swt_merge) <- c("study_id", "glucose_reading", "bp_sys_1", "bp_dia_1",
                      "medication", "medication_reason", "hiv_pos", "bmi", "site")
swt_merge$hiv_pos <- ifelse(swt_merge$hiv_pos == "Positive", 1, 0)

# bushbuckridge 1
bbr_1_merge <- bbr_1 %>%
  filter(!is.na(study_id) & study_id != "") %>%
  select(-c(record_id, data_capturer, datacatureinitial_complete)) %>%
  mutate(site = "Bushbuckridge")

names(bbr_1_merge)[names(bbr_1_merge) == "medicine_use"] <- "medication"
names(bbr_1_merge)[names(bbr_1_merge) == "consent_hiv_test"] <- "consent_hiv"
names(bbr_1_merge)[names(bbr_1_merge) == "collection1_date"] <- "date"

# bushbuckridge 2
bbr_2_merge <- bbr_2 %>%
  select(-c(microbiome2018_complete)) %>%
  mutate(site = "Bushbuckridge")

names(bbr_2_merge) <- c("study_id", "date", "consent_study",
                        "consent_blood_pressure", "bp_sys_1", "bp_dia_1",
                        "bp_sys_2", "bp_dia_2", "bp_sys_3", "bp_dia_3",
                        "height", "weight", "bmi", "antibiotic_use_last",
                        "diarrhea_last", "medication", "consent_hiv","hiv_pos",
                        "consent_glucose", "glucose_reading", "instructions",
                        "instructions_clear", "ashamed_participating",
                        "future_participate", "comments", "site")

# add data for sample in MB17094
MB17094_meta <- voke_meta %>%
  filter(SampleID == "MB17094") %>%
  select(SampleID, bmi_values, hiv_result, glucose_reading) %>%
  mutate(site = "Bushbuckridge")

names(MB17094_meta) <- c("study_id", "bmi", "hiv_pos",
                         "glucose_reading", "site")

pheno_merge <- bind_rows(swt_merge, bbr_1_merge, bbr_2_merge, MB17094_meta)

# add sample names
pheno_merge$sample <- sequenced_samples[match(pheno_merge$study_id,
                                              sequenced_samples$study_id), "sample"]

pheno_merge <- pheno_merge %>%
  filter(study_id %in% sequenced_samples$study_id) %>%
  distinct(study_id, .keep_all = TRUE)

# how many HIV neg individuals excluded? ----
pheno_merge %>%
  group_by(site) %>%
  tally()

pheno_merge %>%
  group_by(hiv_pos, site) %>%
  tally()

# keep HIV negative individuals and remove repeat MB17051 entry ----
pheno_merge_hiv_neg <- pheno_merge %>%
  filter(hiv_pos == 0) %>%
  distinct(study_id, .keep_all = TRUE)

# add sequence id
pheno_merge_hiv_neg <- pheno_merge_hiv_neg %>%
  mutate(sample = sequenced_samples[match(study_id, sequenced_samples$study_id),
                                    "sample"]) %>%
  select(sample, study_id, site, date, hiv_pos, height, weight, bmi,
         glucose_reading, bp_sys_1, bp_dia_1, bp_sys_2, bp_dia_2, bp_sys_3,
         bp_dia_3, starts_with("medication"), ends_with("last"),
         starts_with("consent"), starts_with("instructions"),
         ashamed_participating, future_participate, comments)

pheno_merge_hiv_neg <- pheno_merge_hiv_neg[mixedorder(pheno_merge_hiv_neg$sample), ]

# bmi group ----
# set outlier to NA
# add classification -- lean, overweight, obese
pheno_merge_hiv_neg <- pheno_merge_hiv_neg %>%
  mutate(bmi = ifelse(bmi < 15, NA, bmi),
         bmi_status = ifelse(bmi < 25, "Lean",
                             ifelse(bmi < 30, "Overweight", "Obese")))
pheno_merge_hiv_neg$bmi_status <- factor(pheno_merge_hiv_neg$bmi_status,
                                         levels = c("Lean", "Overweight", "Obese"))

# meta and pheno dfs ----
za_pheno <- pheno_merge_hiv_neg
za_meta <- za_pheno %>%
  select(sample, study_id, site)

# add age ----
# participant ages from AWI-Gen
mb_ages <- read.table(here("input_final/pheno/mbiom_ages.tsv"),
                      sep = "\t", col.names = c("sample1", "sample2",
                                                "awi_date", "awi_age"))

mb_ages_long <- mb_ages %>%
  pivot_longer(cols = c(sample1, sample2))

# awigen to mb link
mb_pheno <- read.csv(here("input_final/pheno/mb_pheno.csv"))

# figure out mapping
za <- merge(za_meta, mb_pheno[, c("mb_id", "fid")], by.x = "study_id",
            by.y = "mb_id", all.x = T)
za <- merge(za, mb_pheno[, c("mb_id", "fid")], by.x = "study_id",
            by.y = "fid", all.x = T)

za_long <- za %>%
  pivot_longer(cols = c(study_id, fid, mb_id)) %>%
  left_join(mb_ages_long, by = "value")

## add in one missing Soweto age -- 1333396
soweto_mb_age_ht_wt <- read.table(here("input_final/pheno/soweto_microbiome_age_ht_wt.txt"),
                                  sep = "\t", header = T)

za_long[which(za_long$value == "1333396"), "awi_age"] <-
  floor(as.numeric(soweto_mb_age_ht_wt[which(soweto_mb_age_ht_wt$MicrobiomeStudy_ID == "1333396"), "AGE"]))

# add dummy date
za_long[which(za_long$value == "1333396"), "awi_date"] <- "2017-10-01"

dates <- za_long %>%
  filter(!is.na(awi_age)) %>%
  select(sample, site, awi_date, awi_age) %>%
  distinct()

dates$awi_date <- as.Date(dates$awi_date, format = "%Y-%m-%d")

# calculate age from date of enrollment
pheno_date <- merge(dates, za_pheno, by = c("sample", "site"), all.x = T)
pheno_date$date <- as.Date(gsub("\\s.+", "", pheno_date$date), format = "%m/%d/%y")

# no date for Soweto - call it 10/1/17
pheno_date$date[is.na(pheno_date$date)] <- as.Date("2017-10-01")

pheno_date <- pheno_date %>%
  mutate(
    age_mb_exact = awi_age + (as.numeric(date - awi_date) / 365),
    age_mb = awi_age + floor(as.numeric(date - awi_date) / 365)
  )

za_pheno <- pheno_date %>%
  mutate(age = age_mb) %>%
  select(-c(awi_date, awi_age, age_mb_exact, age_mb))

# df <- pheno_date %>%
#   select(study_id, site, awi_date, date, awi_age, age_mb_exact, age_mb)
# names(df) <- c("sample_id", "site", "awigen_date", "mb_date", "awigen_age", "mb_age_fraction", "mb_age")
# write.table(df, "/Users/tamburif/Desktop/mb_pilot_ages.txt", sep = "\t", row.names = F, quote = F)

## label samples something helpful for publishing
za_pheno <- za_pheno %>%
  mutate(site_abbrev = ifelse(site == "Bushbuckridge", "BBR", "SWT")) %>%
  arrange(site, sample) %>%
  group_by(site) %>%
  mutate(label = paste(tolower(site), row_number(), sep = "_"),
         label_abbrev = paste(site_abbrev, row_number(), sep = "")) %>%
  ungroup()

labels <- za_pheno %>%
  select(sample, study_id, site, label, label_abbrev)

write.table(labels, here("input_final/za_labels.tsv"), sep = "\t",
            row.names = F, quote = F)

# save pheno data ----
# saveRDS(za_pheno, here("rds/za_pheno.rds"))
# saveRDS(za_meta, here("rds/za_meta.rds"))

save(za_meta, za_pheno, file = here("RData/metadata.RData"))

# setdiff(voke_meta$SampleID, pheno_merge_hiv_neg$study_id)
# setdiff(pheno_merge_hiv_neg$study_id, voke_meta$SampleID)
