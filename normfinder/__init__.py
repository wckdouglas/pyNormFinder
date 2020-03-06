#!/usr/bin/env python

import pandas as pd
import numpy as np
from collections import Counter
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('NormFinder')


def normalize_count(count_mat, return_sf = False):
    '''
    DESeq2 size factor:
        https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106#Sec22
    input:
        pandas dataframe with shape(m,n): m is gene count, n is sample number
        return_sf: return size factors (n)
    output:
        np.ndarray(m,n)
    '''

    cols = count_mat.columns.tolist()
    log_row_geomean = count_mat.transform(np.log)\
            .mean(axis=1, skipna=True)
    finite_mean = np.isfinite(log_row_geomean)
    
    
    def get_size_factor(sample_cnts):
        ratio = np.log(sample_cnts) - log_row_geomean
        ratio = ratio[(finite_mean) & (sample_cnts>0)]
        ratio = np.median(ratio)
        return ratio
    
    sf = np.apply_along_axis(get_size_factor, 0, count_mat)
    sf = np.exp(sf)
    norm_count = count_mat/sf
    logger.info('Performed DESeq2 normalization')
    if return_sf:
        return norm_count, sf
    else:
        return norm_count


class NormFinder:
    '''
        A python implementation for norm finder:
        https://moma.dk/files/r.NormOldStab5.txt


        input:
            - table pandas data frame: 
                first columns is gene name

                Gene X1 X2
            1   U6   10 20
            2   U5   1  4

            - groups: list
                must be same length as shape(table)[1] - 1
                define the groupings of the samples

            - normalized: if the data is already normalized?
                if False, a deseq2 normalization will be applied
            - log: if true, it will perform log2 transform the data


        example:
            nf = NormFinder(count_table, groups=['group1','group2','group1'...],normalized=False)
            stability_df = nf.calc_stability()
        '''

    def __init__(self, table, groups=None, normalized=False, log=True):
        self.table = table
        self.groups = np.array(groups)
        self.group_factor = set(self.groups)
        self.gene_names = self.table.iloc[:, 0] # first column is always gene name
        self.data = self.table.iloc[:, 1:] # the rest are expression data
        self.ngenes = self.table.shape[0]
        self.ngroups = len(self.group_factor)
        self.sample_names = self.data.columns
        if not normalized:
            self.data = normalize_count(self.data)
        
        assert len(self.groups) == self.data.shape[1], "Supplied groups doesn't match the shape of input table"
        self.group_mapper = {s:g for s,g in zip(self.sample_names, self.groups)}

        if log:
            self.log_data = np.log2(self.data + 0.05)
            logger.info('Performed log transformation')
        else:
            self.log_data = self.data 
        self.log_data = self.log_data  \
            .assign(gene = self.gene_names) \
            .pipe(pd.melt, id_vars = 'gene', 
                            var_name='samplename', 
                            value_name='expression')   \
            .assign(expression = lambda d: d.expression.fillna(0))\
            .pipe(lambda d: d[~pd.isnull(d.expression)])  \
            .assign(group = lambda d: d.samplename.map(self.group_mapper)) 
        logger.info('Prepared input table')

    
    def calc_stability(self):
        '''
        a tidy implementation:
        https://rdrr.io/github/dhammarstrom/generefer/src/R/tidy_normfinder.R

        output:
            stability_df: pandas dataframe with columns:
                1. gene: gene name as defined by the input table
                2. stability: NormFinder stability value
                3. gene_mean: mean expression of the gene
        '''
        if self.ngroups > 1:
            logger.info('Using multigroup: %i groups' %self.ngroups)
            stability_df = self.log_data \
                .groupby('samplename', as_index=False)\
                .apply(lambda d: d.assign(sample_mean = lambda d: d.expression.mean())) \
                .groupby(['group'], as_index=False)\
                .apply(lambda d: d.assign(group_mean = lambda d: d.expression.mean())) \
                .groupby(['group','gene'], as_index=False)\
                .apply(lambda d: d.assign(gene_group_mean = lambda d: d.expression.mean())) \
                .groupby('gene', as_index=False)\
                .apply(lambda d: d.assign(gene_mean = lambda d: d.expression.mean())) \
                .assign(global_mean = lambda d: d.expression.mean())\
                .groupby(['gene', 'group'],as_index=False) \
                .apply(lambda d: d.assign(N = lambda d: d.shape[0])\
                                    .assign(residual = lambda d: np.sum((d.expression - d.gene_group_mean - \
                                                                d.sample_mean + d.group_mean)**2\
                                                                /(d.N-1) ))) \
                .groupby('group', as_index=False) \
                .apply(lambda d: d.assign(var_sum_group = np.sum(d.residual.unique()))) \
                .assign(var = lambda d: (d.residual - d.var_sum_group / \
                                    (self.ngenes * (self.ngenes-1))) / (1-2/self.ngenes) ) \
                .groupby(['gene','group'], as_index=False)\
                .agg({'gene_group_mean':'mean',
                    'group_mean':'mean',
                    'gene_mean':'mean',
                    'global_mean':'mean',
                    'N':'mean',
                    'var':'mean'}) \
                .assign(diff = lambda d: d.gene_group_mean - d.group_mean - d.gene_mean  + d.global_mean) \
                .assign(var_n = lambda d: d['var'] / d.N) \
                .assign(sumdd = lambda d: np.sum(d['diff']*d['diff'])) \
                .assign(meanvar = lambda d: d.var_n.mean()) \
                .assign(n = lambda d: (self.ngroups - 1) * (self.ngenes-1) ) \
                .assign(t = lambda d: np.where(d.sumdd/d.n > 0, d.sumdd/d.n, 0 ) )\
                .assign(dn = lambda d: d['diff'] * d['t'] / (d['t'] + d.var_n ) ) \
                .assign(varnew = lambda d: d['var_n'] + d['t'] * d['var_n']/(d['t'] + d['var_n']) ) \
                .assign(qm = lambda d: np.abs(d.dn) + np.sqrt(d.varnew) ) \
                .groupby('gene', as_index=False) \
                .agg({'qm':'mean', 
                      'gene_mean':'mean'}) \
                .rename(columns={'qm':'stability'}) \
                .sort_values('stability') 

        else:
            logger.info('Using single group')
            stability_df = self.log_data \
                .groupby('sample', as_index=False) \
                .apply(lambda d: d.assign(sample_mean = lambda d: d.expression.mean())) \
                .group_by('gene', as_index=False)\
                .apply(lambda d: d.assign(gene_mean = lambda d: d.expression.mean())) \
                .assign(global_mean = lambda d: d.expression.mean())\
                .groupby('gene', as_index=False) \
                .apply(lambda d: d.assign(N = lambda d: d.shape[0])\
                                    .assign(residual = lambda d: np.sum((d.expression - d.gene_mean - \
                                                                d.sample_mean + d.global_mean)**2\
                                                                /(d.N-1) ))) \
                .assign(var_sum = lambda d: d.residual.unique()/sum()) \
                .assign(var = lambda d: (d.residual - var_sum / (self.ngenes * (self.ngenes -1))) /
                                        (1-2/self.ngenes) )\
                .groupby('gene', as_index=False)\
                .agg({'var': lambda x: np.sqrt(x.mean()),
                    'gene_mean': 'mean'}) \
                .rename(columns={'var':'stability'}) \
                .sort_values('stability') 
    
        return stability_df


