# Melissa: Bayesian clustering and imputation of single cell methylomes

New technologies enabling the measurement of DNA methylation at the single cell level are promising to revolutionise our understanding of epigenetic control of gene expression. Yet, intrinsic limitations of the technology result in very sparse coverage of CpG sites (around 5% to 20% coverage), effectively limiting the analysis repertoire to a semi-quantitative level.

__Melissa__ (MEthyLation Inference for Single cell Analysis), is a Bayesian hierarchical method to quantify spatially-varying methylation profiles across genomic regions from single-cell bisulfite sequencing data (scBS-seq). Melissa clusters individual cells based on local methylation patterns, enabling the discovery of epigenetic diversities and commonalities among individual cells. The clustering also acts as an effective regularisation method for imputation of methylation on unassayed CpG sites, enabling transfer of information between individual cells. 

<!--- ![Melissa model overview](analysis/model/figures/melissa.png) -->

<img src="inst/figures/melissa.png" alt="" style="width: 50px;"/> 


The probabilistic graphical representation of the Melissa model is shown below:

![](inst/figures/melissa-model-small.png)


## Installation
To get the latest development version from Github:

```R
# install.packages("devtools")
devtools::install_github("andreaskapou/Melissa", build_vignettes = TRUE)
```

### Melissa dependence

Melissa depends heavily on the [BPRMeth package](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bty129), which is available on 

`Bioconductor`: [http://bioconductor.org/packages/BPRMeth/](http://bioconductor.org/packages/BPRMeth/) and 

`Github`: [https://github.com/andreaskapou/BPRMeth](https://github.com/andreaskapou/BPRMeth).

### Archive repository
There is also an archived version of the Melissa model for reproducing the results presented in the Kapourani and Sanguinetti (2018) biorXiv paper [https://github.com/andreaskapou/Melissa-archive](https://github.com/andreaskapou/Melissa-archive).


## Citation
Kapourani, C.-A. and Sanguinetti, G. (2018). Melissa: Bayesian clustering and imputation of single cell methylomes, bioRxiv.
