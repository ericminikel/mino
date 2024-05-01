options(stringsAsFactors=F)
setwd('~/d/sci/src/mino/')
library(tidyverse)
library(janitor)
library(openxlsx)

# FIGURE 1 stuff
tmt_pep = read.xlsx('data/iqp/1622_deliverables/1622_P1_Peptides.xlsx', startRow=3) %>%
  as_tibble() %>%
  clean_names()
write_tsv(tmt_pep, 'data/tmt/tmt_1622_peptides.tsv')


# FIGURE 3 stuff
sex_meta = tibble(sex=c('F','M'),pch=c(16,17), color=c('#D71234','#1234D7'))
write_tsv(sex_meta, 'data/animals/sex_meta.tsv')

main = read_tsv('data/animals/CMR-2285-full.tsv', col_types = cols()) %>%
  clean_names() %>%
  mutate(animal = as.character(id)) %>%
  rename(sex=gender) %>%
  select(animal, dob, sex)

write_tsv(main, 'data/animals/CMR-2285-main.tsv')
  
ellas = rbind_files('~/d/sci/src/ella/data/export/','.csv') %>%
  mutate(plate = as.integer(substr(file,5,8)))

ellas %>%
  mutate(animal = substr(sample_name,1,8)) %>%
  filter(animal %in% main$animal) %>%
  mutate(date_text = substr(sample_name,10,19)) %>%
  mutate(date = case_when(nchar(date_text) < 10 ~ as.Date(date_text,'%m-%d-%Y'),
                          nchar(date_text) == 10 ~ as.Date(date_text,'%Y-%m-%d'))) %>%
  mutate(nfl = as.numeric(gsub('<','',mean_conc))) %>%
  select(plate, inlet, sample_name, animal, date, gn_rs, nfl) -> cmr2285_ella

write_tsv(cmr2285_ella, 'data/ella/CMR-2285-ella.tsv')




# FIGURE 4 stuff

# list of microglia DEGs for mouse cortex, exported from Cloupe Browser,
# from Mortberg & Gentile 2023 PMID: 37188501 
mg_genes = read_csv('data_received/microglia_genes.csv') %>%
 clean_names()
write_tsv(mg_genes, 'data/tmt/microglia_genes.tsv')

tmt_pep = read.xlsx('data_received/1997_deliverables/1997_P1_peptide.xlsx', startRow=3) %>%
  as_tibble() %>%
  clean_names()
write_tsv(tmt_pep, 'data/tmt/tmt_1997_peptides.tsv')

ellas %>%
  mutate(plate = as.integer(substr(file,5,8))) %>%
  filter(plate==18) %>%
  filter(grepl('^MD', sample_name)) %>%
  mutate(animal = as.integer(gsub('-','',substr(sample_name, 4,5)))) %>%
  mutate(nfl = as.numeric(gsub('<','',mean_conc))) %>%
  select(plate, inlet, sample_name, animal, gn_rs, nfl) -> ella_ef0048dame

write_tsv(ella_ef0048dame, 'data/ella/EF-0048-DA-ME-ella.tsv')




ellas %>%
  mutate(plate = as.integer(substr(file,5,8))) %>%
  filter(plate %in% 11:16) %>%
  filter(grepl('^MP', sample_name)) %>%
  arrange(desc(plate)) %>%
  filter(!duplicated(sample_name)) %>% # use only re-runs where applicable
  mutate(animal = as.integer(gsub('-','',substr(sample_name, 4,5))),
         day = ifelse(grepl('PreTx',sample_name),0,round(as.numeric(gsub('.+-','',sample_name))/24))) %>%
  mutate(nfl = as.numeric(gsub('<','',mean_conc))) %>%
  select(plate, inlet, sample_name, animal, day, gn_rs, nfl) -> ella_ef0047dame
write_tsv(ella_ef0047dame, 'data/ella/EF-0047-DA-ME-ella.tsv')


read_tsv('data/animals/CMR-1982-master.tsv', col_types=cols()) %>%
  mutate(animal = as.character(animal)) %>%
  mutate(comments = case_when(exp_comments == 'fight wound on back' ~ 'fight wound',
                              TRUE ~ '')) %>%
  mutate(tx = case_when(tx=='n/a' ~ 'none',
                        TRUE ~ tx)) %>%
  select(animal, sex, tx, dose, comments) -> main
write_tsv(main, 'data/animals/CMR-1982-main.tsv')

ellas %>%
  mutate(animal=gsub('_.*','',sample_name)) %>%
  mutate(day = case_when(grepl('plasma-pre',sample_name) ~ 0,
                         grepl('plasma$',sample_name) ~ 5)) %>%
  mutate(nfl = as.numeric(gsub('<','',mean_conc))) %>%
  mutate(plate = as.integer(substr(file,5,8))) %>%
  filter(plate %in% 30:31) %>%
  semi_join(main, by='animal') %>%
  select(animal, day, gn_rs, nfl) -> cmr1982ella
write_tsv(cmr1982ella, 'data/ella/CMR-1982-ella.tsv')



ef0047wt_raw = read.xlsx('data_received/EF-0047-DA-ME Final -In Life Data Summary.xlsx',sheet='Dose Administration',startRow=8,
                         colNames = F) %>% 
  as_tibble()
colnames(ef0047wt_raw) = c('day','date','drug','group','animal','route','weight','formulation_admin_g',
                           'conc_mgml','density_gml','dose_administered_mg','dose_administered_mgkg','protocol_specified_mgkg',
                           'variance_percent')
ef0047wt_raw %>%
  fill(day) %>%
  mutate(weight_g = 1000*suppressWarnings(as.numeric(weight))) %>%
  select(day, animal, weight_g) %>%
  filter(!is.na(weight_g)) -> ef0047wt
write_tsv(ef0047wt, 'data/animals/EF-0047-DA-ME-weights.tsv')










main = read_tsv('data/animals/CMR-2535-main-raw.tsv', col_types=cols()) %>%
  clean_names() %>%
  mutate(sex = ifelse(sex,'M','F')) %>% # by default it is read as a T/F logical. convert to M/F character
  mutate(animal = as.character(animal)) %>%
  mutate(tx = gsub(' .*','',tx)) %>%
  mutate(bleed1_vol = as.numeric(substr(bleed1_comments,1,2)),
         bleed2_vol = as.numeric(bleed2_comments)) %>%
  mutate(bleed_vol = case_when(!is.na(bleed1_vol) ~ bleed1_vol,
                               !is.na(bleed2_vol) ~ bleed2_vol,
                               TRUE ~ 80)) %>%
  mutate(bleed_day = case_when(!is.na(bleed1_date) ~ 5,
                         !is.na(bleed2_date) ~ 10)) %>%
  select(animal, bleed_day, tx, sex, bleed_vol) %>%
  mutate(cohort = paste0(tx,'-',bleed_day))
write_tsv(main, 'data/animals/CMR-2535-main.tsv')

weights = read_tsv('data_received/CMR-2535-weights.tsv', col_types=cols()) %>%
  select(-LATEST) %>%
  pivot_longer(cols=-animal) %>%
  mutate(date=as.Date(name), weight=value) %>%
  mutate(day = as.integer(date - min(date) + 1)) %>%
  filter(!is.na(weight)) %>%
  group_by(animal) %>%
  mutate(delta = weight / weight[day==1] - 1) %>%
  select(animal, day, weight, delta)
write_tsv(weights, 'data/animals/CMR-2535-weights.tsv')

ellas %>% 
  mutate(plate = as.integer(substr(file, 5, 8))) %>%
  filter(plate %in% c(65,66,68)) %>%
  mutate(nfl = as.numeric(gsub('<','',mean_conc))) %>%
  mutate(animal = substr(sample_name,1,8)) %>%
  mutate(day = case_when(animal %in% main$animal[main$bleed_day==5] ~ 5,
                         grepl('_mino wash out bleed',sample_name) ~ 24,
                         animal %in% main$animal[main$bleed_day==10] ~ 10)) %>%
  filter(animal %in% main$animal) %>%
  select(animal, day, gn_rs, nfl) -> cmr_2535_ella
write_tsv(cmr_2535_ella, 'data/animals/CMR-2535-ella.tsv')

deg = read_tsv('output/minocycline_saline.tsv', col_types=cols()) %>%
  rename(l2fc = log2fold_change) %>%
  select(-enr)
write_tsv(deg, 'data/rnaseq/mp_deseq2_dge.tsv')         


# FIGURE 5 STUFF
coc_main = read_tsv('data_received/coculture_nfl.tsv', col_types=cols()) %>%
  clean_names() %>%
  fill(plate, group, id) %>%
  mutate(group = case_when(grepl('Co-Mino2[5-7]',group) ~ 'Co-Mino25', # fix someone's Excel drag-to-fill errors
                           TRUE ~ group)) %>%
  filter(group %in% c('Co-CTRL','Co-Mino25')) %>%
  mutate(treatment = case_when(group == 'Co-CTRL' ~ 'control',
                               group == 'Co-Mino25' ~ 'minocycline 25 µM')) %>%
  select(plate, treatment, id, day=timepoint, nfl=nf_l_pg_ml)
write_tsv(coc_main, 'data/cultures/coc.tsv')

