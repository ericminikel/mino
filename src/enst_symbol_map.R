library(tidyverse)
library(janitor)
library(stringr)
gtf = read_tsv('~/d/sci/src/mino/data/rnaseq/ref/gtf_short.txt', col_names = F) %>% clean_names()
gtf$transcript = gsub(';.*','',gtf$x1)
gtf$symbol = gsub('gene_name ','',gsub(';','',str_extract(gtf$x1, 'gene_name.*;')))

gtf %>%
  select(transcript, symbol) %>%
  group_by(transcript, symbol) %>%
  slice(1) -> gtf_uniq

write_tsv(gtf_uniq, '~/d/sci/src/mino/data/rnaseq/ref/gtf_uniq.tsv')
