install.packages("remotes")
remotes::install_github("oganm/ogbox") #for softParser
install.packages("BiocManager")
library(ogbox)
library(BiocManager)
BiocManager::install("GEOquery")
library("GEOquery")
install.packages("matrixStats")
library(matrixStats)

library(tidyverse)
library(devtools)
library(ggplot2)
library(here)
library(stringr)



# read in labonte expr data
labonte_expr = read.csv(file = here("data", "GSE102556", "GSE102556_HumanMDD_fpkmtab.txt.gz"), sep = '\t')

soft_file = here("data", "GSE102556", "GSE102556_family.soft")

labonte_meta = softParser(soft_file)

labonte_meta_pfc = labonte_meta %>% filter(`!Sample_characteristics_ch1 = tissue` == 'Dorsolateral prefrontal cortex (dlPFC; BA8/9)')

gse=getGEO(filename=here("data", "GSE102556", "GSE102556-GPL11154_series_matrix.txt"))
gse_df = gse %>% as.data.frame()

gse_df$title %>% head

sample_nums = str_extract(gse_df$title, '[0-9]+')
gse_df$subject_num = sample_nums
gse

labonte_expr %>% colnames

#new_expr_colnames = colnames(labonte_expr[, 2:282]) %>% str_extract(., 'X[0-9]+')
#colnames(labonte_expr)[2:282] = new_expr_colnames
#labonte_expr

long_region_names = c('Orbitofrontal (OFC; BA11)', 'Dorsolateral prefrontal cortex (dlPFC; BA8/9)', 'Cingulate gyrus 25 (Cg25)', 'Anterior Insula (aINS)', 'Nucleus Accumbens (Nac)', 'Subiculum (Sub)')
short_region_name = c('BA11', 'BA8_9', 'BA25', 'AntIns', 'Nac', 'Subic')

gse_df = gse_df %>%
  mutate(short_region_names = plyr::mapvalues(tissue.ch1, long_region_names, short_region_name))
gse_df = gse_df %>% mutate(age.ch1 = as.numeric(age.ch1))
gse_df = gse_df %>% mutate(expr_names = paste0('X', sample_nums, '.', short_region_names))

labonte_meta = gse_df %>% filter(short_region_names == 'BA8_9')

dlpfc_samples = gse_df %>% filter(short_region_names == 'BA8_9') %>% pull(expr_names)

dlpfc_expr = labonte_expr[c('gene_name', dlpfc_samples)]

labonte_dlpfc_expr = dlpfc_expr

write.csv(labonte_meta, file = here("data", "labonte_dlpfc_meta.csv"))
write.csv(labonte_dlpfc_expr, file = here("data","labonte_dlpfc_expr.csv"))

