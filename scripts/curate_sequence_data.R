library(genefilter)
library(here)
library(metagenomeSeq)
library(tidyverse)

## metadata ----
za_meta <- readRDS("rds/za_meta.rds")
za_pheno <- readRDS("rds/za_pheno.rds")

## public global sequence data ----
# exclude rampelli participants under the age of 18
exclude_rampelli <- c("SRR1929485", "SRR1929574", "SRR1930128", "SRR1930132",
                      "SRR1930134")

pheno_global <- read.table(here("input_final/pheno_global.txt"), sep = "\t",
                           header = T)

pheno_global <- pheno_global %>%
  filter(!(sample %in% exclude_rampelli))

bf_S <- read.table(here("input_final/kraken/genbank_jan2020/burkina_faso/bracken_S_reads.txt"),
                   sep = "\t", header = T, quote = "", comment.char = "")

pheno_bf <- data.frame("sample" = colnames(bf_S), "site" = "Burkina Faso",
                       "site2" = "Burkina Faso")

pheno_global <- rbind(pheno_global, pheno_bf)
pheno_global$site2 <- factor(pheno_global$site2,
                             levels = c("Tanzania", "Madagascar", "Burkina Faso",
                                        "Bushbuckridge", "Soweto", "Sweden",
                                        "United States"))

saveRDS(pheno_global, here("rds/pheno_global.rds"))

### kraken/bracken output ----
for (rank in c("S", "G", "F", "O", "C", "P")){
  
  # za reads
  za_path <- here(paste0("input_final/kraken/genbank_jan2020/za/bracken_",
                    rank, "_reads.txt"))
  za_reads <- read.table(za_path, sep = "\t", check.names = F, quote = "",
                         comment.char = "")
  za_reads <- za_reads[, za_meta$sample]
  
  # global reads
  global_path <- here(paste0("input_final/kraken/genbank_jan2020/global/bracken_",
                  rank, "_reads.txt"))
  global_reads <- read.table(global_path, sep = "\t", check.names = F, quote = "",
                             comment.char = "")
  global_reads <- global_reads[, !grepl("[A-C]\\d+$", names(global_reads))]
  
  # rampelli reads
  rampelli_path <- here(paste0("input_final/kraken/genbank_jan2020/rampelli/bracken_",
                    rank, "_reads.txt"))
  rampelli_reads <- read.table(rampelli_path, sep = "\t", check.names = F, quote = "",
                               comment.char = "")
  
  # burkina faso reads
  bf_path <- here(paste0("input_final/kraken/genbank_jan2020/burkina_faso//bracken_",
                         rank, "_reads.txt"))
  bf_reads <- read.table(bf_path, sep = "\t", check.names = F, quote = "",
                         comment.char = "")
  
  # exclude non-adult samples
  rampelli_reads <- rampelli_reads[, which(!(names(rampelli_reads) %in%
                                               exclude_rampelli))]
  
  # merge datasets
  global_reads <- merge(global_reads, za_reads, by = "row.names",
                        all = T) %>%
    column_to_rownames("Row.names")
  
  global_reads <- merge(global_reads, rampelli_reads, by = "row.names",
                        all = T) %>%
    column_to_rownames("Row.names")
  
  global_reads <- merge(global_reads, bf_reads, by = "row.names",
                        all = T) %>%
    column_to_rownames("Row.names")
  
  global_reads[is.na(global_reads)] <- 0
  
  global_reads <- global_reads[, pheno_global$sample]
  global_reads <- global_reads[!(rownames(global_reads) %in% c("Homo", "Homo sapiens")), ]
  global_reads <- global_reads[rev(order(rowMeans(global_reads))), ]
  
  # filter out rare/low-abundance uninformative features
  global_reads <- global_reads[genefilter(global_reads, pOverA(p = 0.01, A = 0)), ]
  
  # za_reads <- global_reads[, za_meta$sample]
  za_reads <- za_reads[genefilter(za_reads, pOverA(p = 0.01, A = 0)), ]
  za_reads <- za_reads[rev(order(rowMeans(za_reads))), ]
  
  global_rel <- sweep(global_reads, 2, colSums(global_reads), FUN = "/")
  za_rel <- sweep(za_reads, 2, colSums(za_reads), FUN = "/")
  
  global_pseudo <- global_reads + 1
  global_pseudo_rel <- sweep(global_pseudo, 2, colSums(global_pseudo), FUN = "/")
  
  za_pseudo <- za_reads + 1
  za_pseudo_rel <- sweep(za_pseudo, 2, colSums(za_pseudo), FUN = "/")
  
  # save rds
  saveRDS(global_reads, here("rds", paste0("global_", rank, ".rds")))
  saveRDS(global_rel, here("rds", paste0("global_", rank, "_rel.rds")))
  saveRDS(global_pseudo_rel, here("rds", paste0("global_", rank, "_pseudo_rel.rds")))
  
  saveRDS(za_reads, here("rds", paste0("za_", rank, ".rds")))
  saveRDS(za_rel, here("rds", paste0("za_", rank, "_rel.rds")))
  saveRDS(za_pseudo_rel, here("rds", paste0("za_", rank, "_pseudo_rel.rds")))
}

## css normalize
# za
for (rank in c("S", "G")){
  cts <- readRDS(here(paste0("rds/za_", rank, ".rds")))
  mr <- newMRexperiment(cts)
  p <- cumNormStatFast(mr)
  mr_css <- cumNorm(mr, p = p)
  counts_css <- MRcounts(mr_css, norm = T, log = T)
  saveRDS(counts_css, here("rds", paste0("za_", rank, "_css.rds")))
}

# global
for (rank in c("S", "G")){
  cts <- readRDS(here(paste0("rds/global_", rank, ".rds")))
  mr <- newMRexperiment(cts)
  p <- cumNormStatFast(mr)
  mr_css <- cumNorm(mr, p = p)
  counts_css <- MRcounts(mr_css, norm = T, log = T)
  saveRDS(counts_css, here("rds", paste0("global_", rank, "_css.rds")))
}

## sourmash data ----
for (k in c("21", "31", "51")){
  for (suffix in c("", "_track_abund")){
    
    sourmash <- read.csv(here(paste0("input_final/sourmash/compare_k", k,
                                     suffix, ".csv")),
                         header = T)
    
    names(sourmash) <- gsub("X.oak.+_trim_kmers.|_concat.fq.abundtrim", "",
                            names(sourmash))
    
    rownames(sourmash) <- names(sourmash)
    
    # exclude rampelli samples and other extraneous samples
    keep <- which(names(sourmash) %in% pheno_global$sample)
    sourmash <- sourmash[keep, keep]
    
    saveRDS(sourmash, here("rds", paste0("sourmash_k", k, suffix, ".rds")))
  }
}