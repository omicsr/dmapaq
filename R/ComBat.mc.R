#' A multi-processor wrapper for ComBat method.
#'
#' A multi-processor wrapper for ComBat method.
#' ComBat is a method to adjust batch effect where the batch covariate is known.
#'
#' @param dat A data matrix with column for samples and row for probe.
#' @param batch Batch covariate (multiple batches allowed)
#' @param nCores Number of cores will be used for computation
#' @param ... See `ComBat` in `sva` package for extra options
#'
#' @return A data matrix with the same dimension as input data, adjusted for batch effects.
#'     Warning: Values for multimodal distributed CpGs could be over-adjusted.
#'
#' @references Johnson, WE, Rabinovic, A, and Li, C (2007).
#'     *Adjusting batch effects in microarray expression data using Empirical Bayes methods. Biostatistics 8(1):118-127.*
#'
#' @seealso See `ComBat` in `sva` package for extra options
#'
#' @importFrom foreach `%dopar%`
#'
#' @export
ComBat.mc <- function(dat, batch, nCores = 1, ...) {
  #    if(!require("sva")){stop("Can not load sva package")}
  if (nCores > parallel::detectCores()) {
    nCores <- parallel::detectCores()
    cat(
      "Only ", parallel::detectCores(), " Cores avialable, nCores was reset to ",
      parallel::detectCores(), "\n"
    )
  }
  cat("Analysis is running, please wait...!", "\n")

  rname <- rownames(dat)
  nparts1 <- min(ceiling(ncol(dat) / 500), floor(nrow(dat) / 40000))
  parts1 <- rep(1:nparts1, ceiling(nrow(dat) / nparts1))[1:nrow(dat)]
  nCores <- min(10, nCores, ceiling(nrow(dat) / (nparts1 * 20000)))
  dat.o <- NULL
  c1 <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(c1)
  for (i in 1:nparts1) {
    dat1 <- dat[which(parts1 == i), ]
    parts <- rep(1:nCores, ceiling(nrow(dat1) / nCores))[1:nrow(dat1)]
    parts <- sample(parts)
    dat.o1 <- foreach::foreach(s = 1:nCores, .combine = rbind, .export = c("sva::ComBat")) %dopar% {
      s <- s
      idx <- which(parts == s)
      sva::ComBat(dat = dat1[idx, ], batch = batch, ...)
    }
    dat.o <- rbind(dat.o, dat.o1)
  }
  parallel::stopCluster(c1)
  dat.o[rname, ]
}
