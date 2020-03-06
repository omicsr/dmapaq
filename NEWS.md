# dmapaq (development version)

## New features

* `read_idats()` allows to efficiently import idats files mostly 
    using [minfi](https://doi.org/doi:10.18129/B9.bioc.minfi) functions.
* `qc_idats()` allows to compute quality-control of methylation array from Illumina 
    using a [rmarkdown template](inst/rmarkdown/templates/qc_idats/skeleton/skeleton.Rmd).
* `ggheatmap()` allows to compute heatmap with dendrogram on x-axis and y-axis 
    using [ggplot2](https://ggplot2.tidyverse.org/).
* `ComBat.mc()` is a multi-processor wrapper for ComBat method 
    (initially in [ENmix](https://bioconductor.org/packages/ENmix/)).
