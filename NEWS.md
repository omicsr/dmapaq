# dmapaq 0.3.5

* In `DESCRIPTION`,
    + Remove `doParallel`.
* In [rmarkdown templates](inst/rmarkdown/templates),
    + Make sure `Sample_ID` is a unique string and in the first column of the sample sheet.
    + Remove `doParallel`.
    + Order PCA association tile plot according to p-values.
* In `read_idats.R`,
    + Rewrite the build process of the `RGChannelSetExtended` object.

# dmapaq 0.3.4

* In [rmarkdown templates](inst/rmarkdown/templates),
    + Explicit use of `Sample_ID` as the unique IDs for samples.
* In `read_idats.R`,
    + Explicit use of `Sample_ID` as the unique IDs for samples.

# dmapaq 0.3.3

* In [rmarkdown templates](inst/rmarkdown/templates),
    + Update obsolete captions.
    + Perform independent association testing between technical variables and principal components.
    + Fix `ggplot2` issue with legend components defined in `theme_set()`.

# dmapaq 0.3.2

* In [rmarkdown templates](inst/rmarkdown/templates),
    + Remove unnecessary lower condition for outliers.
    + Replace underscore with dash.
    + Fix PCA section.

# dmapaq 0.3.1

* In [rmarkdown templates](inst/rmarkdown/templates),
    + Fix reference for ComBat.
    + Fix missing params in YAML.
* In `read_idats.R`,
    + Fix dataset for population specific CpG/SNPs.
    + Fix missing params in YAML.

# dmapaq 0.3.0

* In [rmarkdown templates](inst/rmarkdown/templates),
    + Reorder code in setup chunk.
    + Rewrite code using `data.table`.
* In `ggheatmap.R`,
    + Remove function `ggheatmap()`.
* In `read_idats.R`,
    + Optimise code.
    + Use `parallel`.
    + Use `data.table`.
    + Add `echo` parameter.
    + Return log messages.

# dmapaq 0.2.0

## Minor improvements and fixes

* In `DESCRIPTION`, 
    + Update all packages listed in `Imports`.
    + Now imports `ggplot2 v3.3.0`.
    + Use `gt` for tables.
* In `README.Rmd`, remove parallel function for ComBat batch correction.
* In [rmarkdown template](inst/rmarkdown/templates/qc_idats/skeleton/skeleton.Rmd),
    + Fix regular expression for reference panel.
    + Fix library availability check for cordbloodcombined reference panel.
    + Fix misplaced parenthesis around calls to `system.file()`.
    + Fix missing probe id in beta matrix output.
    + Fix missing prefix.
    + Fix cell composition section when no tissue is provided.
    + Fix gender discrepancy table when observed gender is missing.
    + Fix call rates table not displaying all thresholds.
    + Remove parallel computation for ComBat batch correction.
    + Replace `requireNamespace()` with `system.file()` in order to check packages availability.
    + No longer show empty gender check table.
* In `R/qc_idats.R`.
    + Import `gt` for tables.
    + Remove import for `kableExtra`.
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
