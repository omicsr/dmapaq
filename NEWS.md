# dmapaq 0.1.3

## Minor improvements and fixes

* Fix missing probe id in beta matrix output from [rmarkdown template](inst/rmarkdown/templates/qc_idats/skeleton/skeleton.Rmd)

# dmapaq 0.1.2

## Minor improvements and fixes

* Now imports `ggplot2 v3.3.0`.
* Update all packages listed in `Imports`.

# dmapaq 0.1.1

## Minor improvements and fixes

* Fix missing prefix in [rmarkdown template](inst/rmarkdown/templates/qc_idats/skeleton/skeleton.Rmd).
* Remove `R/ComBat.mc`.
* Use `parallel` for ComBat in [rmarkdown template](inst/rmarkdown/templates/qc_idats/skeleton/skeleton.Rmd).

# dmapaq 0.1.0

## New features

* `read_idats()` allows to efficiently import idats files mostly 
    using [minfi](https://bioconductor.org/packages/minfi/) functions.
* `qc_idats()` allows to compute quality-control of methylation array from Illumina 
    using a [rmarkdown template](inst/rmarkdown/templates/qc_idats/skeleton/skeleton.Rmd).
* `ggheatmap()` allows to compute heatmap with dendrogram on x-axis and y-axis 
    using [ggplot2](https://ggplot2.tidyverse.org/).
* `ComBat.mc()` is a multi-processor wrapper for ComBat method 
    (initially in [ENmix](https://bioconductor.org/packages/ENmix/)).
