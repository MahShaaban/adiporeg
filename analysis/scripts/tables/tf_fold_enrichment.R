# load required libraries
library(tidyverse)
library(xtable)

# read tf_tf_targets data
tf_targets <- read_csv('data/tf_tf_targets.csv') %>%
  filter(grepl("5' UTR", annotation) | grepl('Promoter', annotation)) %>%
  select(factor, SYMBOL, fold_enrichment, stage) %>%
  mutate(SYMBOL = toupper(SYMBOL)) %>%
  setNames(c('from', 'to', 'weight', 'stage')) %>%
  na.omit()

tf_targets %>%
  group_by(stage, from, to) %>%
  summarise(average = round(mean(weight), 2),
            sd = round(sd(weight), 2)) %>%
  mutate(num = ifelse(is.na(sd), as.character(average), paste(as.character(average), '$\\pm$', as.character(sd)))) %>%
  select(-average, -sd) %>%
  ungroup() %>%
  spread(stage, num) %>%
  mutate(from = ifelse(duplicated(from), '', from)) %>%
  setNames(c('TF', 'Gene', 'Stage 0', 'Stage 1', 'Stage 3')) %>%
  xtable(caption = 'Fold enrichment of transcription factors genes.',
         label = 'tab:tf_fold_enrichment',
         align = 'cllccc') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        add.to.row = list(pos = list(6, 12, 19, 25, 31, 36),
                          command = rep('\\midrule ', 6)),
        caption.placement = 'top',
        sanitize.text.function = identity,
        comment = FALSE,
        file = 'manuscript/tables/tf_fold_enrichment.tex')
