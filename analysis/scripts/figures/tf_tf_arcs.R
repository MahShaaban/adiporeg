# load required libraries
library(tidyverse)
library(reshape2)
library(tidygraph)
library(ggraph)

# load enrichment data
tf_targets <- read_csv('data/tf_tf_targets.csv')

df <- tf_targets %>%
  filter(grepl("5' UTR", annotation) | grepl('Promoter', annotation)) %>%
  select(factor, SYMBOL, fold_enrichment, stage, annotation) %>%
  mutate(SYMBOL = toupper(SYMBOL)) %>%
  setNames(c('from', 'to', 'weight', 'stage', 'annotation')) %>%
  na.omit()

(df %>%
    group_by(from, to, stage) %>%
    summarise(weight = mean(weight)) %>%
    igraph::graph_from_data_frame() %>%
    as_tbl_graph() %>%
    ggraph(layout = 'linear', circular = TRUE) +
    geom_edge_arc(aes(width = weight, color = weight), arrow = arrow()) +
    geom_node_text(aes(label = name)) +
    facet_wrap(~stage) +
    theme_graph()) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/tf_tf_arcs.png',
         width = 25, height = 10, units = 'cm')

