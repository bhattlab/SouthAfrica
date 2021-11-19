# Bushbuckridge and Soweto microbiome project

Analysis scripts for Wits-Stanford shotgun microbiome data.

To generate all figures from the manuscript, please follow these steps:

1. Clone this repo and create a new R project (in RStudio, File > New Project > Existing Directory)
2. Access phenotype and medication data from EGA under dataset ID [EGAD00001006581](https://ega-archive.org/datasets/EGAD00001006581). You can submit an access request to the H3Africa Data Access Committee at https://catalog.h3africa.org
3. Edit `scripts/create_metadata.R` to point to the locations of the EGA phenotype and medication data and run this script to create the required RData file
4. Install all required packages and run each R script in `scripts/figures` and `scripts/figures/supplementary`

     e.g. from the command line:

          find scripts/figures -name "*.R" | xargs -I foo Rscript foo

Required packages:

* DESeq2
* MASS
* Maaslin2
* RColorBrewer
* cowplot
* dplyr
* genefilter
* ggplot2
* ggpubr
* ggrepel
* ggspatial
* gplots
* gtools
* harrietr
* here
* maps
* metagenomeSeq
* reshape2
* rgeos
* rnaturalearth
* rnaturalearthdata
* scales
* sf
* tidyverse
* vegan
