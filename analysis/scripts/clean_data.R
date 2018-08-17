# load required libraries
library(tidyverse)

# load enrichment data
if(!file.exists('data/tf_tf_targets.csv')) {
  runs <- read_csv('data/chip.csv') %>%
  full_join(read_csv('data/chip_fastq.csv')) %>%
  filter(factor %in% c('CEBPA', 'CEBPB', 'CTCF', 'EP300', 'MED1', 'PPARG', 'RXRG')) %>%
  select(factor, stage, run)

  annotated_peaks_files <- list.files('data/annotated_peaks/', pattern = '*.csv', full.names = TRUE)
  names(annotated_peaks_files) <- str_split(annotated_peaks_files, '/|\\.', simplify = TRUE)[, 4]
  enrichment <-  annotated_peaks_files %>%
    map(function(x) {
      read_csv(x) %>%
        filter(SYMBOL  %in% str_to_title(runs$factor))
    }) %>%
    bind_rows(.id = 'run')

  full_join(runs, enrichment) %>%
    write_csv('data/tf_tf_targets.csv')
}
