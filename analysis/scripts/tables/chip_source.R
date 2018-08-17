# load required libraries
library(tidyverse)
library(xtable)

# read data
df <- read_csv('data/chip.csv')

df %>%
  group_by(bibtexkey) %>%
  summarise(gsm = paste(id, collapse = ', '),
            factor = paste(unique(factor), collapse = ', ')) %>%
  mutate(bibtexkey = paste0('\\cite{', bibtexkey, '}')) %>%
  setNames(c('Ref', 'GSM', 'Factors')) %>%
  xtable(caption = 'Sources of the ChIP-Seq data',
         label = 'tab:chip_sources',
         align = 'ccp{.55\\textwidth}p{.35\\textwidth}') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        sanitize.text.function = identity,
        comment = FALSE,
        file = 'manuscript/tables/chip_source.tex')
