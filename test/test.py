%load_ext autoreload
%autoreload 2

import pandas as pd
import numpy as np
from pyNormFinder import NormFinder

def make_data():
    df = pd.read_excel('test_data/CPM_All_Samples.xlsx')   \
        .nlargest(200, '30102-002 - CPM')\
        .drop('Identifier', axis=1)\
        .to_csv('test_data/test_cpm.csv', index=False)      

def make_R(df, groups):
    df.set_index('Name').pipe(lambda d: d+0.05 ).to_csv('test_data/test_R.tsv', sep='\t')
    group_df = {s:g for s,g in zip(df.columns[1:], groups)}
    with open('test_data/test_R.tsv', 'a') as f:
        pd.DataFrame.from_dict({'group':group_df}).transpose().to_csv(f, sep='\t', header=False)



df = pd.read_csv('test_data/test_cpm.csv') \
    .rename(columns=lambda x: x.replace(' ',''))
groups = np.repeat(['KO','CTRL'], 10)
#make_R(df, groups)


R_validate = pd.read_csv('test_data/normfinder_out.csv')
nf = NormFinder(df, groups=groups, normalized=True) 
stab = nf.calc_stability()