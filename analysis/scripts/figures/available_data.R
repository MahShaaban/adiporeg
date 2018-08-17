# load required libraries
library(tidyverse)

# read data
df <- read_csv('data/chip.csv')

(df %>%
  group_by(factor_type, factor, stage) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = factor, y = as.factor(stage), size = n)) +
  geom_point() +
  facet_wrap(~factor_type, ncol = 1, scales = 'free_x', strip.position = 'right') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = 'Differentiation Stage')) %>%
  ggsave(plot = ., 'manuscript/figures/available_data.png')
