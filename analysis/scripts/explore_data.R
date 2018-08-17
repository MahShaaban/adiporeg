# load libraries
library(tidyverse)

# TF with full time points
## extract ids
ind <- read_csv('data/chip.csv') %>%
  filter(factor %in% c('CEBPA', 'CEBPB', 'CTCF', 'EP300', 'MED1', 'PPARG', 'RXRG'),
         stage != 2) %>%
  pull(id)

## get url for fastq files
#read_csv('data/chip_fastq.csv') %>%
#  filter(type == 'SINGLE', id %in% ind) %>%
#  pull(ftp) %>%
#  write('data/tf_full.urls')

