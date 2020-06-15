# read in labonte expr data

library(tidyverse)
library(devtools)
devtools::install_github('oganm/ogbox')
library(ogbox)
library(ggplot2)

labonte_expr = read.csv(file = '~/Downloads/GSE102556_HumanMDD_fpkmtab.txt', sep = '\t')

soft_file = '~/Downloads/GSE102556_family.soft'

labonte_meta = softParser(soft_file)

labonte_meta_pfc = labonte_meta %>% filter(`!Sample_characteristics_ch1 = tissue` == 'Dorsolateral prefrontal cortex (dlPFC; BA8/9)')

source("http://bioconductor.org/biocManager")
biocLite("GEOquery")

install.packages("BiocManager")
library(BiocManager)
BiocManager::install("GEOquery")

library("GEOquery")

gse=getGEO(filename="~/Downloads/GSE102556-GPL11154_series_matrix.txt")
gse_df = gse %>% as.data.frame()

gse_df$title %>% head

library(stringr)

sample_nums = str_extract(gse_df$title, '[0-9]+')
gse_df$subject_num = sample_nums
gse

labonte_expr %>% colnames

new_expr_colnames = colnames(labonte_expr[, 2:282]) %>% str_extract(., 'X[0-9]+')
colnames(labonte_expr)[2:282] = new_expr_colnames
labonte_expr

long_region_names = c('Orbitofrontal (OFC; BA11)', 'Dorsolateral prefrontal cortex (dlPFC; BA8/9)', 'Cingulate gyrus 25 (Cg25)', 'Anterior Insula (aINS)', 'Nucleus Accumbens (Nac)', 'Subiculum (Sub)')
short_region_name = c('BA11', 'BA8_9', 'BA25', 'AntIns', 'Nac', 'Subic')

gse_df = gse_df %>%
  mutate(short_region_names = plyr::mapvalues(tissue.ch1, long_region_names, short_region_name))
gse_df = gse_df %>% mutate(age.ch1 = as.numeric(age.ch1))
gse_df = gse_df %>% mutate(expr_names = paste0('X', sample_nums, '.', short_region_names))

dlpfc_samples = gse_df %>% filter(short_region_names == 'BA8_9') %>% pull(expr_names)


dlpfc_expr = labonte_expr[c('gene_name', dlpfc_samples)]

use_gene_list = c('SST', 'SLC17A7', 'XIST', 'PVALB', 'TAC1', 'MOG', 'GAD1', 'ELFN1')
gene_mat = dlpfc_expr[dlpfc_expr$gene_name %in% use_gene_list, ]
gene_names = gene_mat$gene_name
gene_mat_trans = gene_mat[-1] %>% t() %>% as.data.frame()
colnames(gene_mat_trans) = gene_names

gene_mat_trans = gene_mat_trans %>% tibble::rownames_to_column(var = 'expr_names')

gene_mat_comb = merge(gse_df, gene_mat_trans, by = 'expr_names')

gene_mat_comb %>% ggplot(aes(x = age.ch1, y = ELFN1, color = phenotype.ch1)) + geom_point()

formula = 'SST ~ phenotype.ch1 + age.ch1 + 1 + gender.ch1'
summary(lm(formula, data = gene_mat_comb, na.action = na.omit))
