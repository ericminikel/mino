options(stringsAsFactors=F)
setwd('~/d/sci/src/mino/')
suppressMessages(library(tidyverse))

patient = read_tsv('data/bwh/deid_patient_data.tsv', col_types=cols())

bin_age = function(age) {
  case_when(age < 20 ~ '<20',
            age >= 20 & age < 80 ~ paste0(floor(age/5)*5,'-',floor(age/5)*5+4),
            age >= 80 ~ '≥80')
}

patient %>%
  mutate(age = bin_age(age)) -> patient_binned_age

write_tsv(patient_binned_age, 'data/bwh/deid_patient_data.tsv')
