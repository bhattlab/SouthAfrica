library(dplyr)
library(gtools)
library(here)

# metadata
za_meta <- readRDS(here("rds/za_meta.rds"))

# labels
labels <- read.table(here("input_final/za_labels.tsv"), sep = "\t", header = T)

za_meta$id <- labels[match(za_meta$sample, labels$sample), "id"]

# awigen id
awigen_ids <- read.table(here("input_final/awigen_labels_SRA.txt"), sep = "\t")

## match mb id to awigen id
# mismatched ids
za_meta$study_id[!(za_meta$study_id %in% awigen_ids$V1)]

# fix mismatched ids
mismatched_ids <- za_meta$study_id[!(za_meta$study_id %in% awigen_ids$V1)]

# mismatched ids which don't need to be fixed
mismatched_ids[mismatched_ids %in% awigen_ids$V2]

# fix only mismatched ids
# replace only if in first column (id to be replaced)
za_meta$awigen_id <- ifelse(za_meta$study_id %in% awigen_ids$V1,
                            awigen_ids[match(za_meta$study_id, awigen_ids$V1), "V2"],
                            za_meta$study_id)

za_meta$awigen_id[za_meta$awigen_id == 1333396] <- NA

za_meta <- za_meta[mixedorder(za_meta$id), ]

# get collection dates
pheno_date <- read.table(here("input_final/pheno/pheno_date.txt"),
                         sep = "\t", header = T, quote = "")

sra_data <- za_meta %>%
  mutate(sample_name = id,
         organism = "human gut metagenome",
         host = "Homo sapiens",
         isolation_source = "Human stool",
         collection_date = pheno_date[match(sample, pheno_date$sample), "date"],
         geo_loc_name = ifelse(site == "Bushbuckridge",
                               "South Africa: Bushbuckridge, Mpumalanga",
                               "South Africa: Soweto, Gauteng"),
         lat_lon = ifelse(site == "Bushbuckridge",
                          "24.82 S 31.26 E",
                          "26.25 S 27.85 E")) %>%
  dplyr::select(sample_name, organism, host, isolation_source, collection_date,
         geo_loc_name, lat_lon, awigen_id)

write.table(sra_data, here("SRA_data.txt"), sep = "\t",
            row.names = F, quote = F)

sra_metadata <- za_meta %>%
  mutate(sample_name = id,
         library_ID = sample,
         title = "WGS sequencing of human gut metagenome",
         library_strategy = "WGS",
         library_source = "METAGENOMIC",
         library_selection = "RANDOM",
         library_layout = "PAIRED",
         platform = "ILLUMINA",
         instrument_model = "Illumina HiSeq 4000",
         design_description = "DNA extracted from stool via bead beating, \
         libraries prepared using Illumina Nextera XT, sequenced 150bp PE reads \
         on Illumina HiSeq 4000",
         filetype = "fastq",
         filename = paste0(id, "_R1.fastq"),
         filename2 = paste0(id, "_R2.fastq"))

write.table(sra_metadata, here("SRA_metadata.txt"), sep = "\t", row.names = F, quote = F)