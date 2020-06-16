
library(tidyverse)
library(ggplot2)

library(matrixStats)
library(cowplot)

theme_set(theme_cowplot())

# read expression meta and expr
labonte_meta = read.csv(file = 'data/labonte_dlpfc_meta.csv')
labonte_expr = read.csv(file = 'data/labonte_dlpfc_expr.csv')
labonte_expr = labonte_expr[-1]

#devtools::install_github('oganm/markerGeneProfile', force = T)
library(markerGeneProfile)

# read marker genes
acc_markers = read_csv('marker_genes/ACC_results.csv')
cell_types = acc_markers$adapted_cluster_name %>% unique()

marker_list = lapply(cell_types, function(cell_type){
  print(cell_type)
  return(acc_markers %>% filter(adapted_cluster_name == cell_type) %>% pull(gene) %>% unlist())
})
names(marker_list) = cell_types


# process expression matrix prior to MGP analysis
use_gene_list = labonte_expr %>% pull(gene_name) %>% unlist %>% unique()
gene_mat = labonte_expr[labonte_expr$gene_name %in% use_gene_list, ]
gene_mat = gene_mat %>% distinct(gene_name, .keep_all = T)
gene_names = gene_mat$gene_name
rownames(gene_mat) = rownames(gene_mat) %>% make.names(unique=T)
gene_mat_trans = gene_mat[-1] %>% t() %>% as.data.frame()
colnames(gene_mat_trans) = gene_names

# remove genes with low standard deviations
gene_sds = rowSds(gene_mat[-1] %>% as.matrix())
names(gene_sds) = gene_names
gene_mat = gene_mat[gene_sds > .1, ]
gene_mat_trans = gene_mat_trans %>% tibble::rownames_to_column(var = 'expr_names')

# merge gene expression and meta data frames
gene_mat_comb = merge(labonte_meta, gene_mat_trans, by = 'expr_names')

# plot SST mRNA vs age
p2 = gene_mat_comb %>% ggplot(aes(x = age.ch1, y = SST, color = phenotype.ch1, group = 1)) +  
  geom_smooth(method = "lm", se = F) + geom_point() + 
  ylab('SST mRNA (FPKM)') + xlab('Age (years)')

# run MGP analysis
estimations =  mgpEstimate(exprData=gene_mat,
                           genes=marker_list,
                           geneColName='gene_name',
                           outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                           geneTransform =NULL, # this is the default option for geneTransform
                           groups= NULL, #if there are experimental groups provide them here. if not desired set to NULL
                           seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                           removeMinority = TRUE) 

mgp_estimates = as.data.frame(estimations$estimates) %>% tibble::rownames_to_column(var = 'expr_names')

mgp_df = merge(labonte_meta, mgp_estimates)

# plot SST cell type proportion (MGP) vs age
p1 = mgp_df %>% ggplot(aes(x = age.ch1, y = SST, color = phenotype.ch1, group = 1)) + 
  geom_smooth(method = "lm", se = F) + geom_point() + 
  ylab('SST cell type proportion (MGP)') + xlab('Age (years)')

# plot subplots together
plot_grid(p2, p1, nrow = 1)

