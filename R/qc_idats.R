#' Compute quality-control of methylation array from Illumina using a rmarkdown template.
#'
#' @param csv_file A `character`. The path to a CSV file, *i.e.*, a sample sheet describing the data.
#' @param data_directory A `character`. The path to the data directory.
#' @param array A `character`. The array name, *i.e.*, `"EPIC"` or `"450k"`.
#' @param annotation A `character`. The name and version of the annotation package to be used.
#' @param cohort_name A `character`. The name of the studied cohort / population.
#' @param output_file A `character`. The name of the html file produced.
#' @param output_directory A `character`. The path to the output directory.
#' @param filter_snps A `logical`. Should the probes in which the probed CpG falls near a SNP
#'     (according to ([Zhou et al., 2016](https://www.doi.org/10.1093/nar/gkw967))) be removed?
#'     Default is `TRUE`.
#' @param filter_non_cpg A `logical`. Should the non-cg probes be removed?
#' @param filter_xy A `logical`. Should the probes from X and Y chromosomes be removed?
#'     Default is `TRUE`.
#' @param filter_multihit A `logical`. Should the probes which align to multiple locations
#'     (according to [Nordlund et al., 2013](https://www.doi.org/10.1186/gb-2013-14-9-r105)) be removed?
#'     Default is `TRUE`.
#' @param filter_beads A `logical`. Should the probes with a beadcount less than 3 be removed?
#'     Default is `TRUE`.
#' @param population A `character`. Name of the ethnicity population to be used.
#'     Default is `NULL` for none.
#' @param bead_cutoff A `numeric`. The threshold for beadcount.
#'     Default is `0.05`.
#' @param detection_pvalues A `numeric`. The threshold for the detection pvalues above which,
#'     values are considered as missing. Default is `0.01`.
#' @param filter_callrate A `logical`. Should the data be filtered based on call rate metric?
#'     Default is `TRUE`.
#' @param callrate_samples A `numeric`. The call rate threshold for samples, under which samples are excluded.
#'     Default is `0.99`.
#' @param callrate_probes A `numeric`. The call rate threshold for probes, under which probes are excluded.
#'     Default is `1`.
#' @param gender_threshold A `numeric`. The threshold value to discrimate gender based
#'     on sexual chromosomes methylation.
#'    Default is `-2`.
#' @param gender_colname A `character`. The name of the column containing the gender
#'     in the file provided in `csv_file`.
#'     Default is `NULL`.
#' @param norm_background A `character`. Optional method to estimate background normal distribution parameters.
#'     This must be one of `"oob"`, `"est"` or `"neg"`.
#'     Default is `"oob"`.
#' @param norm_dye A `character`. Dye bias correction, "mean": correction based on averaged red/green ratio;
#'     or `"RELIC"`: correction with RELIC method;
#'     or `"none"`: no dye bias correction.
#'     Default is `"RELIC"`.
#' @param norm_quantile A `character`. The quantile normalisation to be used.
#'     This should be one of `"quantile1"`, `"quantile2"`, or `"quantile3"`.
#'     Default is `"quantile1"`.
#' @param cell_tissue A `character`. The cell tissue to be used for cell composition estimation,
#'     using a reference panel (*i.e.*, `"blood"` or `"cordblood"`) or a mathematical deconvolution.
#'     Default is `NULL`.
#' @param pca A `logical`. Whether or not a PCA should be performed on the dataset.
#'     Default is `TRUE`.
#' @param pca_vars A `vector(character)`. Variables to be used with factorial planes.
#'     Default is `c("Sample_Plate", "Sentrix_ID")`.
#' @param max_labels A `numeric`. The maximum number of labels to show on plots. Default is `15`
#' @param title A `character`. The report's title. Default is `paste(array, "Array Quality-Control")`.
#' @param author_name A `character`. The author's name to be printed in the report.
#'     Default is `Unknown`.
#' @param author_affiliation A `character`. The affiliation to be printed in the report.
#'     Default is `NULL`.
#' @param author_email A `character`. The email to be printed in the report.
#'     Default is `NULL`.
#' @param encoding A `character`. The encoding to be used for the html report.
#'     Default is `"UTF-8"`.
#' @param ... Parameters to pass to `rmarkdown::render()`.
#'
#' @return NULL
#'
#' @export
qc_idats <- function(
  csv_file,
  data_directory,
  array = "EPIC",
  annotation = "ilm10b4.hg19",
  cohort_name = "COHORT",
  output_file = paste(cohort_name, array, "QC.html", sep = "_"),
  output_directory = ".",
  filter_snps = TRUE,
  filter_non_cpg = TRUE,
  filter_xy = TRUE,
  filter_multihit = TRUE,
  filter_beads = TRUE,
  population = NULL,
  bead_cutoff = 0.05,
  detection_pvalues = 0.01,
  filter_callrate = TRUE,
  callrate_samples = 0.99,
  callrate_probes = 1,
  gender_threshold = -2,
  gender_colname = NULL,
  norm_background = "oob",
  norm_dye = "RELIC",
  norm_quantile = "quantile1",
  cell_tissue = NULL,
  pca = TRUE,
  pca_vars = c("Sample_Plate", "Sentrix_ID"),
  max_labels = 15,
  title = paste(array, "Array Quality-Control"),
  author_name = "Unknown",
  author_affiliation = NULL,
  author_email = NULL,
  encoding = "UTF-8",
  ...
) {
  message_prefix <- "[dmapaq] "

  message(message_prefix, "Quality-Control started ...")

  invisible(file.copy(
    from = system.file("rmarkdown", "templates", "qc_idats", "skeleton.Rmd", package = "dmapaq"),
    to = file.path(tempdir(), "qc_idats.Rmd"),
    overwrite = TRUE
  ))

  if (!all(sapply(
    X = strsplit(gsub("\\([^()]+\\)|\n| ", "", utils::packageDescription("dmapaq", fields = "Suggests")), ",")[[1]],
    FUN = function(x) nchar(system.file(package = x)) > 0
  ))) {
    stop(message_prefix, "Packages listed in 'Suggests' field are needed to run this function.")
  }

  rmarkdown::render(
    input = file.path(tempdir(), "qc_idats.Rmd"),
    output_file = output_file,
    output_dir = output_directory,
    encoding = encoding,
    params = list(
      csv_file = csv_file,
      data_directory = data_directory,
      array = array,
      annotation = annotation,
      filter_snps = filter_snps,
      filter_non_cpg = filter_non_cpg,
      filter_xy = filter_xy,
      filter_multihit = filter_multihit,
      filter_beads = filter_beads,
      population = population,
      bead_cutoff = bead_cutoff,
      detection_pvalues = detection_pvalues,
      filter_callrate = filter_callrate,
      callrate_samples = callrate_samples,
      callrate_probes = callrate_probes,
      gender_threshold = gender_threshold,
      gender_colname = gender_colname,
      norm_background = norm_background,
      norm_dye = norm_dye,
      norm_quantile = norm_quantile,
      cell_tissue = cell_tissue,
      pca = pca,
      pca_vars = pca_vars,
      max_labels = max_labels,
      title = title,
      author_name = author_name,
      author_affiliation = author_affiliation,
      author_email = author_email,
      output_directory = output_directory
    ),
    ...
  )

  message(message_prefix, "Quality-Control ended.")

  invisible()
}
