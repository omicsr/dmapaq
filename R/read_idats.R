#' Efficiently import idats files mostly using minfi functions.
#'
#' @param directory A `character`. Location of IDAT files, default is the current working directory.
#' @param csv_file A `character`. Path to the sample sheet (csv files) or
#'     name of the sample sheet in `directory`.
#' @param meth_value_type A `character`. Indicates whether you prefer m-values (`"M"`)
#'     or beta-values (`"B"`). Default is `"B"`.
#' @param array_name A `character`. Choose microarray type, eiyther `"450K"` or `"EPIC"`.
#'     Default is `"EPIC"`.
#' @param annotation_version A `character`. Version of the annotation package that should be used.
#'     Default is `"ilm10b4.hg19"` for the `"EPIC"` array
#' @param n_cores An `integer`. The number of cores to use,
#'     i.e., at most how many child processes will be run simultaneously.
#' @param rgSet A `RGChannelSet` object.
#' @param echo A `logical`. Should messages be displayed?
#'
#' @inheritParams qc_idats
#'
#' @import data.table
#'
#' @return A `list`.
#' @export
read_idats <- function(
  directory = getwd(),
  csv_file = "csv$",
  meth_value_type = "B",
  filter_beads = TRUE,
  bead_cutoff = 0.05,
  filter_non_cpg = TRUE,
  filter_snps = TRUE,
  population = NULL,
  filter_multihit = TRUE,
  filter_xy = TRUE,
  detection_pvalues = 0.01,
  filter_callrate = TRUE,
  callrate_samples = 0.99,
  callrate_probes = 1,
  norm_background = "oob",
  norm_dye = "RELIC",
  norm_quantile = "quantile1",
  array_name = c("EPIC", "450k"),
  annotation_version = c("ilm10b4.hg19", "ilmn12.hg19"),
  n_cores = 1,
  rgSet = NULL,
  echo = FALSE
) {

  if (echo) {
    message(
      "===============================\n",
      "[dmapaq] Reading IDAT files ...\n",
      "===============================",
      appendLF = TRUE
    )
  }

  array_name <- array_name[1]
  annotation_version <- annotation_version[1]

  stopifnot(suppressPackageStartupMessages(
    nchar(system.file(package = "minfi")) > 0 &
      switch(
        EXPR = array_name,
        "450k" = nchar(system.file(package = "IlluminaHumanMethylation450kmanifest")) > 0,
        "EPIC" = nchar(system.file(package = "IlluminaHumanMethylationEPICmanifest")) > 0
      )
  ))

  if (is.null(rgSet) | !inherits(rgSet, "RGChannelSet")) {
    sample_sheet <- read_sample_sheet(directory = directory, csv_file = csv_file, echo = echo)
    rgSet <- read_metharray_exp(sample_sheet = sample_sheet, n_cores = n_cores)
  }
  rgSet@annotation <- switch(
    EXPR = array_name,
    "450k" = c(array = "IlluminaHumanMethylation450k", annotation = annotation_version),
    "EPIC" = c(array = "IlluminaHumanMethylationEPIC", annotation = annotation_version)
  )

  data_detP <- minfi::detectionP(rgSet)
  data_detP[is.na(data_detP)] <- 1

  if (filter_callrate) {
    good_detection <- data_detP < detection_pvalues

    call_rate_samples <- colSums(good_detection) / nrow(good_detection)
    bad_samples <- names(which(call_rate_samples < callrate_samples))

    good_detection <- good_detection[, setdiff(colnames(good_detection), bad_samples)]

    call_rate_cpg <- rowSums(good_detection) / ncol(good_detection)
    bad_cpgs <- names(which(call_rate_cpg < callrate_probes))
  } else {
    bad_samples <- NULL
    bad_cpgs <- NULL
  }

  mset_raw <- mset <- minfi::preprocessRaw(rgSet)

  if (filter_beads) {
    bc <- get_beadcount(rgSet)
    bc2 <- bc[rowSums(is.na(bc)) < bead_cutoff * (ncol(bc)), ]
    mset_f2 <- mset[minfi::featureNames(mset) %in% rownames(bc2), ]
    n_beads_discarded <- format(dim(mset)[1] - dim(mset_f2)[1], big.mark = ",", digits = 0)
    mset <- mset_f2
  }

  if (filter_non_cpg) {
    mset_f2 <- minfi::dropMethylationLoci(mset, dropCH = TRUE)
    n_non_cpg_discarded <- format(dim(mset)[1] - dim(mset_f2)[1], big.mark = ",", digits = 0)
    mset <- mset_f2
  }

  if (filter_snps) {
    ref_population <- c(
      "AFR", "EAS", "EUR",
      "SAS", "AMR", "GWD", "YRI", "TSI", "IBS",
      "CHS", "PUR", "JPT", "GIH", "CHB", "STU",
      "ITU", "LWK", "KHV", "FIN", "ESN", "CEU",
      "PJL", "ACB", "CLM", "CDX", "GBR", "BEB",
      "PEL", "MSL", "MXL", "ASW"
    )

    if (is.null(population) || !(population %in% ref_population)) {
      manifest_hg19 <- switch(
        EXPR = array_name,
        "450k" = get(utils::data("hm450.manifest.hg19", package = "ChAMPdata")),
        "EPIC" = get(utils::data("EPIC.manifest.hg19", package = "ChAMPdata"))
      )
      which_population <- which(manifest_hg19$MASK_general)
    } else {
      manifest_hg19 <- switch(
        EXPR = array_name,
        "450k" = get(utils::data("hm450.manifest.pop.hg19", package = "ChAMPdata")),
        "EPIC" = get(utils::data("EPIC.manifest.pop.hg19", package = "ChAMPdata"))
      )
      which_population <- which(manifest_hg19[, paste("MASK_general", population, sep = "_")])
    }
    maskname <- rownames(manifest_hg19)[which_population]
    mset_f2 <- mset[!minfi::featureNames(mset) %in% maskname, ]
    n_snps_discarded <- format(dim(mset)[1] - dim(mset_f2)[1], big.mark = ",", digits = 0)
    mset <- mset_f2
  }

  if (filter_multihit) {
    multi_hit <- get(utils::data("multi.hit", package = "ChAMPdata"))
    mset_f2 <- mset[!minfi::featureNames(mset) %in% multi_hit$TargetID, ]
    n_multihit_discarded <- format(dim(mset)[1] - dim(mset_f2)[1], big.mark = ",", digits = 0)
    mset <- mset_f2
  }

  if (filter_xy) {
    switch(
      EXPR = array_name,
      "450k" = utils::data("probe.features", package = "ChAMPdata"),
      "EPIC" = utils::data("probe.features.epic", package = "ChAMPdata")
    )
    probe_features <- get("probe.features")
    autosomes <- probe_features[!probe_features$CHR %in% c("X", "Y"), ]
    mset_f2 <- mset[minfi::featureNames(mset) %in% rownames(autosomes), ]
    n_xy_discarded <- format(dim(mset)[1] - dim(mset_f2)[1], big.mark = ",", digits = 0)
    mset <- mset_f2
  }

  mset <- ENmix::preprocessENmix(
    rgSet = rgSet,
    bgParaEst = norm_background,
    dyeCorr = norm_dye,
    QCinfo = NULL,
    exQCsample = FALSE,
    exQCcpg = FALSE,
    exSample = bad_samples,
    exCpG = unique(c(bad_cpgs, setdiff(minfi::featureNames(mset_raw), minfi::featureNames(mset)))),
    nCores = n_cores
  )
  mset <- ENmix::norm.quantile(mdat = mset, method = norm_quantile)

  methylation_matrix <-switch(meth_value_type,
    "B" = minfi::getBeta(mset, "Illumina"),
    "M" = minfi::getM(mset),
    stop('Methylation value type not defined. Only "B" or "M" are available.')
  )

  if (min(methylation_matrix, na.rm = TRUE) <= 0) {
    methylation_matrix[methylation_matrix <= 0] <- min(methylation_matrix[methylation_matrix > 0], na.rm = TRUE)
  }
  if (max(methylation_matrix, na.rm = TRUE) >= 1) {
    methylation_matrix[methylation_matrix >= 1] <- max(methylation_matrix[methylation_matrix < 1], na.rm = TRUE)
  }

  colnames(methylation_matrix) <- minfi::pData(mset)[["Sample_ID"]]
  mset@metadata[[meth_value_type]] <- methylation_matrix
  tmp_phenotypes <- data.table::as.data.table(minfi::pData(rgSet))
  tmp_phenotypes[, "Sample_ID" := lapply(.SD, as.character), .SDcols = "Sample_ID"]
  tmp_means <- colMeans(data_detP)[tmp_phenotypes[["Sample_ID"]]]
  tmp_phenotypes[, "mean_detection_pvalue" := tmp_means]
  tmp_callrate <- (colSums(data_detP < detection_pvalues) / nrow(data_detP))[tmp_phenotypes[["Sample_ID"]]]
  tmp_phenotypes[, "call_rate" := tmp_callrate]
  mset@metadata[["phenotypes"]] <- tmp_phenotypes

  log_msg <- character(0)
  if (filter_callrate) {
    log_msg <- c(log_msg, paste0(
      "Filtering samples with call rate below ",
      paste(format(callrate_samples * 100, digits = 1, nsmall = 1), "%"), ":\n",
      "  - ", format(length(bad_samples), big.mark = ",", digits = 0), " samples were discarded."
    ))
    log_msg <- c(log_msg, paste0(
      "Filtering probes with call rate below ",
      paste(format(callrate_probes * 100, digits = 1, nsmall = 1), "%"), ":\n",
      "  - ", format(length(bad_cpgs), big.mark = ",", digits = 0), " probes were discarded."
    ))
  }
  if (filter_beads) {
    log_msg <- c(log_msg, paste0(
      "Filtering probes with a beadcount lower than three in at least ",
      paste(format(bead_cutoff * 100, digits = 1, nsmall = 1), "%"), " of samples:\n",
      "  - ", n_beads_discarded, " probes were discarded."
    ))
  }
  if (filter_non_cpg) {
    log_msg <- c(log_msg, paste0(
      "Filtering non-cg probes:\n",
      "  - ", n_non_cpg_discarded, " probes were discarded."
    ))
  }
  if (filter_snps) {
    log_msg <- c(log_msg, paste0(
      "Filtering probes with SNPs (Zhou et al., 2016; doi:10.1093/nar/gkw967):\n",
      "  - ", n_snps_discarded, " probes were discarded."
    ))
  }
  if (filter_multihit) {
    log_msg <- c(log_msg, paste0(
      "Filtering probes that align to multiple locations ",
      "(Nordlund et al., 2013; doi:10.1186/gb-2013-14-9-r105):\n",
      "  - ", n_multihit_discarded, " probes were discarded."
    ))
  }
  if (filter_xy) {
    log_msg <- c(log_msg, paste0(
      "Filtering probes on the X or Y chromosome:\n",
      "  - ", n_xy_discarded, " probes were discarded."
    ))
  }
  log_msg <- c(log_msg,
    "Zeros have been replaced with smallest value over zero.",
    "Ones have been replaced with largest value below one.",
    paste0(
      "Data contains:\n",
      "  - ", format(dim(methylation_matrix)[1], big.mark = ",", digits = 0), " probes.\n",
      "  - ", format(dim(methylation_matrix)[2], big.mark = ",", digits = 0), " samples.\n",
      "  - ", format(sum(is.na(methylation_matrix)), big.mark = ",", digits = 0), " missing values."
    )
  )
  if (echo) {
    message(paste(log_msg, collapse = "\n"), appendLF = TRUE)
  }

  list(mset = mset, rgset = rgSet, log = log_msg)
}


#' read_sample_sheet
#'
#' @inheritParams read_idats
#' @param ignore.case A `logical`. A logical value.
#'     If `TRUE`, the directory path is prepended to the file names to give a relative file path.
#'     If `FALSE`, the file names (rather than paths) are returned.
#' @param recursive A `logical`. Should the listing recurse into directories?
#' @param full.names A `logical`. Should pattern-matching be case-insensitive?
#'
#' @keywords internal
read_sample_sheet <- function(
  directory = ".",
  csv_file = "csv$",
  ignore.case = TRUE,
  recursive = TRUE,
  full.names = TRUE,
  echo = FALSE
) {
  if (file.exists(suppressWarnings(normalizePath(csv_file)))) {
    list_files <- normalizePath(csv_file)
  } else {
    list_files <- list.files(
      path = directory,
      pattern = csv_file,
      full.names = full.names,
      ignore.case = ignore.case,
      recursive = recursive
    )
    if (length(list_files)>1) {
      warning("[dmapaq] More than one CSV file have been found!")
      list_files <- list.files[1]
      if (echo) message("[dmapaq] File '", list.files, "' will be used.", appendLF = TRUE)
    }
  }

  data_header <- grep("^\\[DATA\\]", readLines(list_files), ignore.case = TRUE)
  if (length(data_header) == 0) data_header <- 0
  col_names <- colnames(data.table::fread(file = list_files, skip = data_header, nrows = 1))
  default_cols <- c("Sample_ID", "Sentrix_ID", "Sentrix_Position")

  cols_missing <- default_cols[!default_cols %in% col_names]
  if (length(cols_missing) != 0) {
    stop(
      "[dmapaq] Sample Sheet must contains the following missing columns:\n",
      "  - ", paste(cols_missing, collapse = "\n  - ")
    )
  }

  sample_sheet <- data.table::fread(file = list_files, skip = data_header)
  data.table::setnames(
    x = sample_sheet,
    old = c("Sample_ID", "Slide", "Array", "Sample_Plate", "Sample_Well"),
    new = c("Sample_ID", "Sentrix_ID", "Sentrix_Position", "Sample_Plate", "Sample_Well"),
    skip_absent = TRUE
  )

  basenames <- sub("_Grn\\.idat.*", "", sapply(
    X = paste0(sample_sheet[["Sentrix_ID"]], "_", sample_sheet[["Sentrix_Position"]], "_Grn.idat"),
    FUN = grep,
    x = list.files(path = directory, recursive = recursive, full.names = TRUE),
    value = TRUE,
    USE.NAMES = FALSE
  ), ignore.case = TRUE)
  sample_sheet[, "Basename" := basenames]

  sample_sheet
}

#' read_metharray
#'
#' @inheritParams read_idats
#' @param files The `basenames` or `filenames` of the IDAT files.
#'     `basenames` are the filename without the ending `_Grn.idat` or `_Red.idat`.
#'     `filenames` are filenames including `_Grn.idat` or `_Red.idat`.
#'
#' @keywords internal
read_metharray <- function(files, n_cores) {
  basenames <- unique(sub("_Grn\\.idat.*|_Red\\.idat.*", "", files))
  for (ichannel in c("Grn", "Red")) {
    i_files <- paste0(basenames, "_", ichannel, ".idat")
    names(i_files) <- basename(basenames)
    i_files_exists <- file.exists(i_files)
    if (!all(i_files_exists)) {
      i_filesgz_exists <- file.exists(paste0(i_files, ".gz"))
      if (!all(i_filesgz_exists)) {
        stop(
           "The following specified files do not exist:\n",
          "  - ", paste(i_files[!i_files_exists], collapse = "\n  - ")
        )
      }
      i_files <- paste0(i_files, ".gz")
    }
    suppressWarnings({
      i_idats <- parallel::mclapply(X = i_files, mc.preschedule = FALSE, mc.cores = n_cores, FUN = illuminaio::readIDAT)
    })
    assign(x = gsub("(^.).*", "\\1_idats", ichannel), value = i_idats)
  }

  common_addresses <- as.character(Reduce("intersect", lapply(
    X = get("G_idats"),
    FUN = function(x) rownames(x$Quants)
  )))

  out <- minfi::RGChannelSetExtended(
    Red = do.call("cbind", lapply(
      X = get("R_idats"),
      y = common_addresses,
      FUN = function(x, y) x$Quants[y, "Mean"]
    )),
    Green = do.call("cbind", lapply(
      X = get("G_idats"),
      y = common_addresses,
      FUN = function(x, y) x$Quants[y, "Mean"]
    )),
    RedSD = do.call("cbind", lapply(
      X = get("R_idats"),
      y = common_addresses,
      FUN = function(x, y) x$Quants[common_addresses, "SD"]
    )),
    GreenSD = do.call("cbind", lapply(
      X = get("G_idats"),
      y = common_addresses,
      FUN = function(x, y) x$Quants[common_addresses, "SD"]
    )),
    NBeads = do.call("cbind", lapply(
      X = get("G_idats"),
      y = common_addresses,
      FUN = function(x, y) x$Quants[common_addresses, "NBeads"]
    ))
  )

  rownames(out) <- common_addresses

  out
}

#' read_metharray_exp
#'
#' @inheritParams read_idats
#' @inheritParams read_metharray
#' @inheritParams read_sample_sheet
#'
#' @keywords internal
read_metharray_exp <- function(
  directory = NULL,
  sample_sheet = NULL,
  ignore.case = TRUE,
  recursive = TRUE,
  full.names = TRUE,
  n_cores = 1
) {
  if (is.null(sample_sheet)) {
    if (is.null(directory)) directory <- "."
    green_files <- list.files(
      path = directory,
      pattern = "_Grn.idat.*",
      recursive = recursive,
      ignore.case = ignore.case,
      full.names = full.names
    )
    red_files <- list.files(
      path = directory,
      pattern = "_Red.idat.*",
      recursive = recursive,
      ignore.case = ignore.case,
      full.names = full.names
    )

    if (length(green_files) == 0 || length(red_files) == 0) {
      stop("IDAT files must be provided.")
    }

    common_files <- intersect(sub("_Grn.idat.*", "", green_files), sub("_Red.idat.*", "", red_files))

    if (length(common_files) == 0) {
      stop('"Grn" and "Red" idats files must be provided.')
    }

    common_files_green <- paste0(common_files, "_Grn.idat")
    if (!setequal(common_files_green, green_files)) {
      warning(
        "The following files only exists for the green channel:\n",
        "  - ", paste(setdiff(green_files, common_files_green), collapse = "\n  - ")
      )
    }

    common_files_red <- paste0(common_files, "_Red.idat")
    if (!setequal(common_files_red, red_files)) {
       warning(
        "The following files only exists for the red channel:\n",
        "  - ", paste(setdiff(red_files, common_files_red), collapse = "\n  - ")
      )
    }

    rgSet <- read_metharray(common_files, n_cores)
  } else {
    if (all(!grepl("Basename", names(sample_sheet)))) {
      stop('"Basename" must be provided as a column of "sample_sheet".')
    }

    if (is.null(directory)) {
      files <- sample_sheet[["Basename"]]
    } else {
      files <- file.path(directory, basename(sample_sheet[["Basename"]]))
    }

    rgSet <- read_metharray(files, n_cores)

    pD <- methods::as(sample_sheet, "DataFrame")
    pD[["filenames"]] <- files
    rownames(pD) <- colnames(rgSet)
    rgSet@colData <- pD
  }

  minfi::sampleNames(rgSet) <- rgSet[["Sample_ID"]]

  rgSet
}

#' get_beadcount
#'
#' @param x `rgSet`
#'
#' @keywords internal
get_beadcount <- function(x) {
  nb <- minfi::getNBeads(x)
  typeI <- minfi::getProbeInfo(x, type = "I")
  typeII <- minfi::getProbeInfo(x, type = "II")
  locus_names <- minfi::getManifestInfo(x, "locusNames")

  bc_temp <- matrix(
    data = NA_real_,
    ncol = ncol(x),
    nrow = length(locus_names),
    dimnames = list(locus_names, minfi::sampleNames(x))
  )
  bc_temp[typeII$Name, ] <- nb[typeII$AddressA, ]
  bc_temp[typeI$Name, ] <- nb[typeI$AddressB, ]
  bc_temp[typeI$Name, ] <- nb[typeI$AddressA, ]
  bc_temp[which(nb[typeI$AddressA, ] < 3 | nb[typeI$AddressB, ] < 3)] <- NA

  data.frame(bc_temp)
}
