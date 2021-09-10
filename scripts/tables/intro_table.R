# Code for introductory table describing gut microbiome literature through 2020

library(here)
library(tidyverse)

studies <- read.table(here("input_final/microbiome_studies.txt"), sep = "\t",
                   header = T, quote = "", skip = 1)

studies <- studies %>%
  mutate(Study.population = tolower(Study.population),
         includes_adults = grepl("adult", Study.population),
         Sequencing.methodology = tolower(Sequencing.methodology),
         used_metagenomics = grepl("metagenomics", Sequencing.methodology))

# countries
non_african <- c("Spain", "Peru", "")
cs <- paste(studies$Country.of.enrollment, collapse = ",")
cs <- gsub(", ", ",", cs)
n_cs <- stringr::str_split(cs, ",")[[1]]
length(unique(n_cs[!(n_cs %in% non_african)]))

# study population
studies %>%
  group_by(includes_adults) %>%
  tally()

# sequencing method
studies %>%
  group_by(used_metagenomics) %>%
  tally()

# metagenomic studies in adults
studies %>%
  group_by(includes_adults, used_metagenomics) %>%
  tally()

studies %>%
  filter(includes_adults, used_metagenomics) %>%
  select(Title, First.author)
