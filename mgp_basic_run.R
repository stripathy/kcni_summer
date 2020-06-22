library(here)
library(tidyverse)
library(ggplot2)

library(matrixStats)
library(cowplot)
library(broom)


theme_set(theme_cowplot())

# install edgeR 
BiocManager::install("edgeR")
library(edgeR)


#devtools::install_github('oganm/markerGeneProfile', force = T)
library(markerGeneProfile)

# read marker genes
acc_markers = read_csv(here('marker_genes','ACC_results.csv'))
cell_types = acc_markers$cluster %>% unique()

marker_list = lapply(cell_types, function(cell_type){
  print(cell_type)
  return(acc_markers %>% filter(cluster == cell_type) %>% pull(gene) %>% unlist())
})
names(marker_list) = cell_types

# convert labonte_expr to cpms

labonte_cpm = edgeR::cpm(labonte_expr, log = F, prior.count = 0)

# process expression matrix prior to MGP analysis
use_gene_list = labonte_cpm %>% rownames() %>% unlist %>% unique()
gene_mat = labonte_cpm[use_gene_list, ]

# calculate genes with very low standard deviations and remove them from the expression matrix
gene_sds = rowSds(gene_mat %>% as.matrix(), na.rm = T) 
gene_mat = gene_mat[gene_sds > .1, ]

gene_mat = gene_mat %>% as.data.frame() %>% tibble::rownames_to_column(var = 'gene_name')

gene_mat = gene_mat %>% distinct(gene_name, .keep_all = T)
gene_names = gene_mat$gene_name
rownames(gene_mat) = rownames(gene_mat) %>% make.names(unique=T)
gene_mat_trans = gene_mat[-1] %>% t() %>% as.data.frame()
colnames(gene_mat_trans) = gene_names

# remove genes with low standard deviations
# gene_sds = rowSds(gene_mat[-1] %>% as.matrix())
# names(gene_sds) = gene_names
# gene_mat = gene_mat[gene_sds > .1, ]
gene_mat_trans = gene_mat_trans %>% tibble::rownames_to_column(var = 'geo_accession')

# merge gene expression and meta data frames
gene_mat_comb = merge(labonte_meta, gene_mat_trans, by = 'geo_accession')


# run MGP analysis
estimations =  mgpEstimate(exprData=gene_mat,
                           genes=marker_list,
                           geneColName='gene_name',
                           outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                           geneTransform =NULL, # this is the default option for geneTransform
                           groups= NULL, #if there are experimental groups provide them here. if not desired set to NULL
                           seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                           removeMinority = TRUE) 

mgp_estimates = as.data.frame(estimations$estimates) %>% tibble::rownames_to_column(var = 'geo_accession')

mgp_df = merge(labonte_meta, mgp_estimates)

# plot SST mRNA vs age
sst_mrna_vs_age_fig = gene_mat_comb %>% ggplot(aes(x = age, y = SST, color = phenotype, group = 1)) +  
  geom_smooth(method = "lm", se = F) + geom_point() + 
  ylab('SST mRNA (CPM)') + xlab('Age (years)')

# plot SST cell type proportion (MGP) vs age
sst_prop_vs_age_fig = mgp_df %>% ggplot(aes(x = age, y = Inh_SST, color = phenotype, group = 1)) + 
  geom_smooth(method = "lm", se = F) + geom_point() + 
  ylab('SST cell type proportion (MGP)') + xlab('Age (years)')

# plot subplots together
plot_grid(sst_mrna_vs_age_fig, sst_prop_vs_age_fig, nrow = 1)

## 

# plot MOG mRNA and Oligo MGP vs age
mog_mrna_vs_age_fig = gene_mat_comb %>% ggplot(aes(x = age, y = MOG, color = phenotype, group = 1)) +  
  geom_smooth(method = "lm", se = F) + geom_point() + 
  ylab('MOG mRNA (CPM)') + xlab('Age (years)')

oligo_prop_vs_age_fig = mgp_df %>% ggplot(aes(x = age, y = Oligo, color = phenotype, group = 1)) + 
  geom_smooth(method = "lm", se = F) + geom_point() + 
  ylab('Oligo cell type proportion (MGP)') + xlab('Age (years)')

plot_grid(mog_mrna_vs_age_fig, oligo_prop_vs_age_fig, nrow = 1)



# merge MGPs with metadata matrix
labonte_meta_plus_mgps = merge(labonte_meta, mgp_df)

# fit a linear model per cell type proportion against all covariates used in DESeq2 modelling
mod_df_list = lapply(cell_types, function(cell_type_name){
  curr_formula = paste0('scale(', cell_type_name, ') ~ gender + scale(ph) + scale(rin) + scale(pmi) + scale(age) + phenotype')
  curr_mod = lm(curr_formula, data = labonte_meta_plus_mgps)
  
  
  mod_df = tidy(curr_mod)
  mod_df$cell_type = cell_type_name
  return(mod_df)
  
}) %>% bind_rows()
# names(mod_df_list = cell_types[1:12])

mod_df_list$padj = p.adjust(mod_df_list$p.value, method = 'BH')

replace_terms = c('Intercept', 'gender:Male', 'pH', 'RIN', 'PMI', 'age', 'disorder:MDD')
mod_df_list$term = plyr::mapvalues(mod_df_list$term, from = c(mod_df_list$term %>% unique), to = replace_terms)

# print data frame for just MDD beta coefficients
mod_df_list %>% filter(term == 'disorder:MDD')

# print data frame for just age beta coefficients
mod_df_list %>% filter(term == 'age')

# beta coeffs per cell type for phenotype and age effects
beta_plot = mod_df_list %>% filter(term %in% c('disorder:MDD', 'age')) %>% 
  ggplot(aes(x = cell_type, y = estimate)) + geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) + 
  ylab('Std. Beta coeff.') + xlab('Cell type proportions') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  facet_wrap(~term)

beta_plot

