# pyNormFinder #

This is a python reimplementation for [NormFinder](https://moma.dk/files/r.NormOldStab5.txt). The implementation of NormFinder follows a *tidy* version of [NormFinder](https://rdrr.io/github/dhammarstrom/generefer/src/R/tidy_normfinder.R).


## Usage ##

```
from pyNormFinder import NormFinder

# input table #
count_table
```

|    | Name           |   001 - CPM |   002 - CPM |   003 - CPM |   004 - CPM |
|---:|:---------------|------------:|------------:|------------:|------------:|
|  0 | hsa-miR-16-5p  |    242712   |    285619   |    457227   |    363398   |
|  1 | hsa-miR-486-5p |    193649   |    238127   |    316513   |    376722   |
|  2 | hsa-let-7b-5p  |     83115.6 |     84112.4 |    116224   |    126185   |


```
groups = ['control','control','KO','KO']
nf = NormFinder(count_table, groups=groups,normalized=False, log=True)
stability_df = nf.calc_stability()
stability_df.head()
```

|     | gene           |   stability |   gene_mean |
|----:|:---------------|------------:|------------:|
| 101 | hsa-miR-26b-5p |    0.17623  |    12.2995  |
| 129 | hsa-miR-361-3p |    0.185044 |     9.42474 |
|  27 | hsa-miR-128-3p |    0.187065 |     9.50192 |


## Credit ##

Please cite the [original paper](https://www.ncbi.nlm.nih.gov/pubmed/15289330) if you plan to use this package in your publication.

- Andersen C.L., Ledet-Jensen J., Ørntoft T.. Normalization of real-time quantitative RT-PCR data: a model based variance estimation approach to identify genes suited for normalization - applied to bladder- and colon-cancer data-sets. **Cancer Research**. 2004 (64): 5245-5250 