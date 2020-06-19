library(here)
library(tidyverse)
library(ggplot2)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

#DESeq2 analyses
labonte_meta = read_csv(file = here('data','labonte_dlpfc_meta.csv'))
colnames(labonte_meta)
labonte_meta %<>% select(geo_accession, expr_names, age = age.ch1, gender = gender.ch1, pmi = pmi.ch1, rin = rin.ch1, ph = ph.ch1, phenotype = phenotype.ch1)

#count data is needed for DEseq2, we obtained this from biojupies
labonte_expr = read_tsv( here('data',"GSE102556","GSE102556-expression.counts.ARCHS4.txt.zip"))
labonte_expr %<>% select(c("gene_symbol", labonte_meta$geo_accession))

#ensure right order
labonte_meta %<>% mutate(geo_accession = factor(geo_accession, levels = colnames(labonte_expr %>% select(-gene_symbol))))
labonte_meta %<>% arrange(geo_accession)
#quick check for order
sum(labonte_meta$geo_accession == colnames(labonte_expr %>% select(-gene_symbol)))

#convert from tibble to data.frame by adding row names
labonte_expr %<>% column_to_rownames("gene_symbol") %>% as.data.frame()

dds <- DESeqDataSetFromMatrix(countData = labonte_expr,
                              colData = labonte_meta,
                              design= ~ gender + phenotype + ph + rin + pmi + age)

#prefilter data
keep <- rowSums(counts(dds)) >= 10
print(paste("Keeping", sum(keep), "of", nrow(labonte_expr), "genes"))
dds <- dds[keep,]

#compute differential expression - may be slow
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
#get the result for case/control status
res <- results(dds, name="phenotype_MDD_vs_CTRL")
res %<>% as_tibble(rownames = "gene_symbol")
#print results ordered by p-value
res %>% arrange(pvalue)

#examine p-value histogram
hist(res$pvalue)

#examine the top gene
#melt data
melted_expression <- labonte_expr %>% as_tibble(rownames = "gene_symbol") %>% melt(value.name = "count") %>% as_tibble()
melted_expression %<>% dplyr::rename(geo_accession = variable)
#select gene
melted_expression %<>% filter(gene_symbol == "RPL23AP4")
#merge metadata
melted_expression <- inner_join(melted_expression, labonte_meta)
melted_expression %>% arrange(-count) #lots of zeroes
melted_expression %>% group_by(phenotype) %>% summarize(mean_count = mean(count), max_count = max(count))

#Student task option - think of and apply more stringent prefiltering

##########################################################################################
# Test polygenic score associations: correlate risk score with expression of each gene
##########################################################################################

risk_scores <- read_csv(here("data", "day1_PRS", "risk_scores.csv"))
labonte_meta <- inner_join(labonte_meta, risk_scores)

#add in a term for day one PRS to the design
dds <- DESeqDataSetFromMatrix(countData = labonte_expr,
                              colData = labonte_meta,
                              design= ~ gender + phenotype + ph + rin + pmi + age + day1_PRS)
#prefilter as before
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds, name="day1_PRS")
res %<>% as_tibble(rownames = "gene_symbol")
res %>% arrange(pvalue)

#SST is of interest lets look at this gene in depth