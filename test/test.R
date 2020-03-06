#!/usr/bin/env Rscript

source('./normfinder.R')

norm <- Normfinder('test_data/test_R.tsv', ctVal=F)
norm$Ordered %>% tibble::rownames_to_column('miRNA') %>% write_csv('test_data/normfinder_out.csv')