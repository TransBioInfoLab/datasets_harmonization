#' Run BMIQ Normalization on DNA Methylation Data
#'
#' Performs Beta Mixture Quantile (BMIQ) normalization on DNA methylation beta values
#' to correct for probe type bias between Infinium Type I and Type II probes. The
#' function applies normalization to each sample in parallel.
#'
#' @param assay_data A numeric matrix of DNA methylation beta values with probes
#'   as rows and samples as columns. Values should be between 0 and 1.
#' @param probe_df A data frame containing probe annotations. Must include a
#'   column named \code{type12} indicating probe type ("I" for Type I probes,
#'   "II" for Type II probes).
#' @param seed Integer. Random seed for reproducibility of the BMIQ algorithm.
#'   Default is 42.
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'   Default is 16.
#'
#' @return A numeric matrix of normalized beta values with the same dimensions
#'   as the input \code{assay_data}. Row and column names are preserved from
#'   the input matrix.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts sample and probe names from the input matrix
#'   \item Converts probe type annotations to numeric format (1 for Type I, 2 for Type II)
#'   \item Handles edge case: If only one sample is provided, creates a duplicate
#'     to allow BMIQ to run (BMIQ may require multiple samples for some operations)
#'   \item Initializes parallel processing backend using \code{start_parallel()}
#'   \item Applies BMIQ normalization to each sample in parallel using \code{plyr::aaply()}
#'   \item Each sample is normalized using \code{BMIQ_adjust()} with the same seed
#'   \item Shuts down parallel processing backend
#'   \item Transposes result and removes duplicate sample if one was added
#'   \item Restores original row and column names
#' }
#'
#' BMIQ normalization adjusts Type II probe distributions to match Type I probe
#' distributions, reducing technical bias in 450K and EPIC array data.
#'
#' @note
#' This function requires the following:
#' \itemize{
#'   \item \code{BMIQ_adjust()} function to be available in the environment
#'   \item \code{start_parallel()} and \code{stop_parallel()} functions for parallel processing
#'   \item The \code{plyr} package for parallel array operations
#' }
#'
run_bmiq <- function(assay_data, probe_df, seed = 42, cores = 16) {
  samples <- colnames(assay_data)
  probes <- row.names(assay_data)
  type12 <- ifelse(probe_df$type12 == "I", 1, 2)
  if (length(samples) == 1) {
    extra_data <- assay_data
    colnames(extra_data) <- paste0("Extra_", colnames(assay_data))
    assay_data <- cbind(assay_data, extra_data)
  }
  
  parallel_res <- start_parallel(TRUE, cores)
  do_parallel <- parallel_res$parallel
  cluster <- parallel_res$cluster
  
  suppressMessages({
    assay_data <- plyr::aaply(
      assay_data, 2,
      function(x){
        set.seed(seed)
        norm_ls <- BMIQ_adjust(
          beta.v = x,
          design.v = type12,
          plots = FALSE,
          pri = FALSE,
          nfit = 50000
        )
        return (norm_ls$nbeta)
      },.progress = "time", .parallel = do_parallel
    )
  })
  
  stop_parallel(do_parallel, cluster)
  
  assay_data <- t(assay_data)
  if (length(samples) == 1) {
    assay_data <- assay_data[,1,drop = FALSE]
  }
  colnames(assay_data) <- samples
  row.names(assay_data) <- probes
  
  assay_data
}

#' Multi-Step BMIQ Normalization with Error Recovery
#'
#' Performs robust BMIQ normalization by attempting multiple random seeds for
#' samples that fail normalization. This function iteratively runs BMIQ with
#' different seeds until all samples are successfully normalized or all seeds
#' are exhausted.
#'
#' @param beta_data A numeric matrix of DNA methylation beta values with probes
#'   as rows and samples as columns. Values should be between 0 and 1.
#' @param probe_df A data frame containing probe annotations. Must include a
#'   column named \code{type12} indicating probe type ("I" for Type I probes,
#'   "II" for Type II probes).
#' @param seeds Numeric vector of random seeds to try sequentially. Default is
#'   \code{c(1, 2, 4, 8, 16, 32)}. Seeds are tried in order until all samples
#'   are successfully normalized.
#' @param cores Integer. Number of CPU cores to use for parallel processing
#'   within each BMIQ run. Default is 16.
#'
#' @return A numeric matrix of normalized beta values with the same dimensions
#'   as the input \code{beta_data}. Row and column names are preserved from
#'   the input matrix. All samples will be normalized unless all seed attempts
#'   fail for a particular sample.
#'
#' @details
#' Sometimes BMIQ fails on a sample for a given random seed due to numerical
#' instabilities or convergence issues. To avoid issues with this, this function
#' runs BMIQ with multiple random seeds for samples that cause errors, until
#' all samples are normalized.
#'
#' The algorithm works as follows:
#' \enumerate{
#'   \item Initialize result matrix and list of target samples (all samples initially)
#'   \item For each seed in the provided seed vector:
#'     \itemize{
#'       \item Run BMIQ normalization on remaining target samples using current seed
#'       \item Identify samples that still contain missing values (normalization failures)
#'       \item Update result matrix with successfully normalized samples
#'       \item Update target list to only include failed samples
#'       \item Print progress information (timestamp, sample counts)
#'     }
#'   \item Continue until all samples are normalized or seeds are exhausted
#' }
#'
#' This approach is particularly useful for large datasets where occasional
#' normalization failures can occur due to extreme data distributions or
#' numerical edge cases.
#'
#' @section Progress Reporting:
#' The function prints progress information for each seed iteration:
#' \itemize{
#'   \item Current seed being processed
#'   \item Start and end timestamps
#'   \item Number of samples being processed
#'   \item Number of samples still missing after the current iteration
#' }
#'
#' @note
#' This function requires:
#' \itemize{
#'   \item \code{run_bmiq()} function to be available in the environment
#'   \item The seed vector should contain diverse values to maximize chances
#'     of successful normalization across all samples
#'   \item Processing time scales linearly with the number of seed attempts needed
#' }
multistep_bmiq <- function(
    beta_data, probe_df, seeds = c(1, 2, 4, 8, 16, 32), cores = 16
) {
  result_data <- beta_data
  
  targets <- colnames(beta_data)
  
  for (seed in seeds) {
    if (length(targets) > 0) {
      print(sprintf("Start Seed: %d", seed))
      print(Sys.time())
      normed <- run_bmiq(
        assay_data = beta_data[,targets, drop = FALSE],
        probe_df = probe_df,
        seed = seed,
        cores = cores
      )
      
      missing_count <- colSums(is.na(normed))
      missing <- sum(missing_count > 0)
      
      print(sprintf("Samples: %d, Missing: %d", length(targets), missing))
      print(Sys.time())
      
      for (target in targets) {
        result_data[,target] <- normed[,target]
      }
      
      targets <- names(missing_count)[missing_count > 0]
    }
  }
  
  result_data
}