# dmapaq 0.2.0

## Minor improvements and fixes

* In `DESCRIPTION`, 
  + Update all packages listed in `Imports`.
  + Now imports `ggplot2 v3.3.0`.
* In `README.Rmd`, remove parallel function for ComBat batch correction.
* In [rmarkdown template](inst/rmarkdown/templates/qc_idats/skeleton/skeleton.Rmd),
    + Fix regular expression for reference panel.
    + Fix library availability check for cordbloodcombined reference panel.
    + Fix misplaced parenthesis around calls to `system.file()`.
    + Fix missing probe id in beta matrix output.
    + Fix missing prefix.
    + Remove parallel computation for ComBat batch correction.
    + Replace `requireNamespace()` with `system.file()` in order to check packages availability.
* Remove `R/ComBat.mc`.

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
