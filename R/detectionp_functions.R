#' Calculate Detection P-Values for DNA Methylation Data
#'
#' Computes detection p-values for each probe in a DNA methylation dataset
#' using the minfi package. Detection p-values assess the confidence that a
#' probe's signal is significantly above background noise. The calculation
#' is performed in parallel across all samples.
#'
#' @param dnam_data A SummarizedExperiment object containing DNA methylation data
#'   with beta values in the "dnam" assay slot. The object must have probes as
#'   rows and samples as columns.
#'
#' @return A numeric matrix of detection p-values with samples as rows and probes
#'   as columns. Row names are sample identifiers from the input data. Lower
#'   p-values indicate higher confidence in probe detection.
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Extracts the DNA methylation assay matrix from the SummarizedExperiment
#'   \item Initializes parallel processing with 8 cores
#'   \item For each sample (column), computes detection p-values using
#'     \code{minfi::detectionP()} with type "mu+u" (based on mean and unmethylated channels)
#'   \item Combines results into a matrix with samples as rows and probes as columns
#'   \item Stops parallel processing
#'   \item Assigns sample names as row names
#' }
#'
#' Detection p-values are calculated by comparing probe intensities to negative
#' control probes. Values closer to 0 indicate stronger signal detection, while
#' values closer to 1 indicate poor detection (signal not distinguishable from background).
detectionp_values <- function(dnam_data) {
  assay_data <- SummarizedExperiment::assays(dnam_data)$dnam
  
  parallel <- TRUE
  parallel <- start_parallel(parallel, cores = 8)
  
  detp_mat <- plyr::ldply(
    1:ncol(assay_data),
    .fun = function(idx) {
      data <- assay_data[,idx]
      res <- minfi::detectionP(data, type = "mu+u")
      t(res)
    }, .parallel = parallel
  )
  stop_parallel(parallel)
  
  row.names(detp_mat) <- colnames(assay_data)
  
  detp_mat
}

#' Filter Probes Based on Detection P-Value Failures
#'
#' Filters a detection p-value matrix to retain only probes that have at least
#' one sample failing the detection threshold. This is useful for identifying
#' problematic probes that need attention or removal from downstream analyses.
#'
#' @param detp_mat A numeric matrix of detection p-values with samples as rows
#'   and probes as columns, typically produced by \code{detectionp_values()}.
#' @param threshold Numeric. The detection p-value threshold for determining
#'   failure. Probes with p-values >= threshold are considered failed detections.
#'   Default is 0.05. Common values range from 0.01 (stringent) to 0.05 (moderate).
#'
#' @return A numeric matrix containing only probes that have at least one sample
#'   with a detection p-value >= threshold. The matrix maintains the same row
#'   structure (samples) as the input but may have fewer columns (probes).
#'
#' @details
#' The filtering process:
#' \enumerate{
#'   \item Creates a logical matrix where TRUE indicates failed detection
#'     (p-value >= threshold)
#'   \item Calculates the total number of failures per probe (column sums)
#'   \item Retains only probes with at least one failure (totals > 0)
#' }
#'
#' This function is useful for quality control workflows where you want to:
#' \itemize{
#'   \item Identify which probes have detection issues in any samples
#'   \item Focus downstream QC efforts on problematic probes
#'   \item Generate reports of probe-level detection failures
#' }
#'
#' Probes that pass detection in all samples (totals == 0) are removed from
#' the output, as they don't require further attention.
detectionp_filter_regular <- function(detp_mat, threshold = 0.05) {
  detp_res <- detp_mat >= threshold
  totals <- colSums(detp_res)
  detp_mat <- detp_mat[,totals > 0]
  
  detp_mat
}

#' Extract Failed Detection Measurements from Detection P-Value Matrix
#'
#' Identifies and extracts all sample-probe combinations that failed the detection
#' p-value threshold, returning a long-format data frame with sample identifiers,
#' probe identifiers, and the corresponding p-values.
#'
#' @param detp_mat A numeric matrix of detection p-values with samples as rows
#'   and probes as columns, typically produced by \code{detectionp_values()}.
#' @param threshold Numeric. The detection p-value threshold for determining
#'   failure. Measurements with p-values >= threshold are considered failures.
#'   Default is 0.05.
#'
#' @return A data frame with three columns:
#'   \describe{
#'     \item{sample}{Character. Sample identifier for the failed measurement}
#'     \item{cpg}{Character. CpG probe identifier for the failed measurement}
#'     \item{detp}{Numeric. The detection p-value for this sample-probe combination}
#'   }
#'   Each row represents one failed measurement. If no measurements fail the
#'   threshold, returns an empty data frame with the same structure.
#'
#' @details
#' The extraction process:
#' \enumerate{
#'   \item Creates a logical matrix where TRUE indicates failed detection
#'     (p-value >= threshold)
#'   \item Identifies row-column indices of all TRUE values using \code{which()}
#'     with \code{arr.ind = TRUE}
#'   \item Extracts the actual p-values for these positions
#'   \item Maps row indices to sample identifiers and column indices to probe identifiers
#'   \item Joins all information into a single data frame
#' }
#'
#' This function is essential for:
#' \itemize{
#'   \item Creating lists of failed measurements for imputation or removal
#'   \item Quality control reporting and visualization
#'   \item Identifying patterns of failures across samples or probes
#'   \item Documenting which measurements were problematic in your dataset
#' }
#'
#' The returned data frame is in long format (one row per failure), making it
#' suitable for analysis, filtering, and integration with other QC workflows.
detectionp_get_failures <- function(detp_mat, threshold = 0.05) {
  detp_res <- detp_mat >= threshold
  
  target_df <- which(detp_res, arr.ind = TRUE) %>%
    as.data.frame()
  
  inds <- cbind(target_df$row, target_df$col)
  target_df$detp <- detp_mat[inds]
  
  cpg_df <- data.frame(
    col = 1:ncol(detp_res),
    cpg = colnames(detp_res)
  )
  
  sample_df <- data.frame(
    row = 1:nrow(detp_res),
    sample = row.names(detp_res)
  )
  
  target_df <- target_df %>%
    dplyr::left_join(cpg_df, by = "col") %>%
    dplyr::left_join(sample_df, by = "row") %>%
    dplyr::select("sample", "cpg", "detp")
  
  target_df
}

#' Calculate Detection P-Value Failures for a Specific Plate
#'
#' Performs detection p-value analysis on samples from a specific plate or batch,
#' combining calculation, filtering, and failure extraction into a single workflow.
#' This is useful for plate-specific quality control in large studies with multiple
#' processing batches.
#'
#' @param dnam_data A SummarizedExperiment object containing DNA methylation data
#'   with beta values in the "dnam" assay slot and plate/batch identifiers in colData.
#' @param plate Character or numeric. The identifier of the plate to analyze. This
#'   value should match entries in the plate column of colData.
#' @param threshold Numeric. The detection p-value threshold for determining
#'   failure. Default is 0.05.
#' @param plate_col Character. The name of the column in colData that contains
#'   plate identifiers. Default is "Sample_Plate".
#'
#' @return A data frame with three columns:
#'   \describe{
#'     \item{sample}{Character. Sample identifier for the failed measurement}
#'     \item{cpg}{Character. CpG probe identifier for the failed measurement}
#'     \item{detp}{Numeric. The detection p-value for this sample-probe combination}
#'   }
#'   Only includes failures from samples belonging to the specified plate.
#'
#' @details
#' This function combines several steps into a single workflow:
#' \enumerate{
#'   \item Extracts plate identifiers from colData using the specified column name
#'   \item Subsets the data to include only samples from the specified plate
#'   \item Calculates detection p-values using \code{detectionp_values()}
#'   \item Filters to probes with failures using \code{detectionp_filter_regular()}
#'   \item Extracts failure details using \code{detectionp_get_failures()}
#' }
#'
#' This plate-specific approach is valuable for:
#' \itemize{
#'   \item Quality control in multi-batch studies
#'   \item Identifying plate-specific technical issues
#'   \item Parallel processing of large datasets by plate
#'   \item Generating plate-level QC reports
#'   \item Deciding whether to exclude entire plates from analysis
#' }
#'
#' The function allows efficient QC assessment when samples are processed across
#' multiple plates, enabling identification of systematic issues specific to
#' particular processing batches.
detectionp_plate <- function(
    dnam_data, plate, threshold = 0.05, plate_col = "Sample_Plate"
) {
  plates <- dnam_data %>%
    SummarizedExperiment::colData() %>%
    as.data.frame() %>%
    dplyr::pull(as.name(plate_col))
  
  dnam_data <- dnam_data[,plates == plate]
  detp_mat <- detectionp_values(dnam_data)
  detp_mat <- detectionp_filter_regular(detp_mat, threshold = threshold)
  detp_df <- detectionp_get_failures(detp_mat, threshold = threshold)
  
  detp_df
}
