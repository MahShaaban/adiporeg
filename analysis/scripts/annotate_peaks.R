# load required libraries
library(tidyverse)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

# create new directory
if(!dir.exists('data/annotated_peaks')) {
  dir.create('data/annotated_peaks')
}

# load peaks
peak_fils <- list.files('data/bed', pattern = '*.xls', full.names = TRUE)
names(peak_fils) <- str_split(peak_fils, '/|\\_', simplify = TRUE)[, 3]

# create file path for annotated peaks
annotated_peaks_files <- paste0('data/annotated_peaks/', names(peak_fils), '.csv')

# read, annotate and write to files
map2(annotated_peaks_files, peak_fils, function(x, y) {
  if(!file.exists(x)) {
    df <- read.delim(y, skip = 23, sep = '\t')
    gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
    annotatePeak(gr,
                 TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                 annoDb = 'org.Mm.eg.db',
                 level = 'gene') %>%
      as.data.frame() %>%
      write_csv(path = x)
    print(paste(x, 'was created.'))
  }
})
