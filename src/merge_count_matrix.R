# gs://fc-3effa716-6fbc-4071-a224-96c78e5129c7/31e9d9d6-17cd-4bf0-896a-7e936d1f4fee/quantify/01043c44-0c4b-48f2-958c-0d8ba933a7d5/call-salmon_paired_reads/cacheCopy/RP-2586_MP-6-384_left-cx_v1_RNA_OnPrem.quant.sf.gz

library(AnVIL)
library(tidyverse)
library(janitor)

# gsutil_cp('gs://fc-3effa716-6fbc-4071-a224-96c78e5129c7/gtf_uniq.tsv', 'R/output/')

files = gsutil_ls('gs://fc-3effa716-6fbc-4071-a224-96c78e5129c7/31e9d9d6-17cd-4bf0-896a-7e936d1f4fee/quantify/*', recursive=T)

quant_files = files[grepl('sf\\.gz',files)]

for (quant_file in quant_files) {
  gsutil_cp(quant_file, 'R/data/')
}

if(exists('quants')) {
  rm(quants)
}
for (quant_file in list.files('R/data/', full.names = T)) {
  temp = read_tsv(quant_file) %>% clean_names()
  sample_id_long = gsub('R/data//','',gsub('_v1.*','',quant_file))
  sample_id_short = gsub('-$','',substr(sample_id_long, 9, 13))
  temp$id_long = sample_id_long
  temp$sample_id = sample_id_short
  if(exists('quants')) {
    quants = rbind(quants, temp)
  } else {
    quants = temp
  }
}

gtf_uniq = read_tsv('R/output/gtf_uniq.tsv')
quants$transcript = gsub('\\..*','',quants$name)
quants$symbol = gtf_uniq$symbol[match(quants$transcript, gtf_uniq$transcript)]
mean(is.na(quants$symbol))  # 0.56%
quants$symbol[is.na(quants$symbol)] = quants$transcript[is.na(quants$symbol)]

quants %>%
  group_by(symbol, sample_id) %>%
  summarize(.groups='keep', tpmsum = sum(tpm)) %>%
  pivot_wider(names_from=sample_id, values_from=tpmsum) -> qmat

write_tsv(qmat, 'R/output/qmat.tsv')
