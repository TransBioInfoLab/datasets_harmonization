#' Order Expression Matrix by Row Standard Deviations
#'
#' Orders the rows of an expression matrix by their standard deviations in
#' descending order, placing the most variable features (e.g., probes) at the top.
#' This is useful for selecting highly variable features for dimensionality
#' reduction analyses.
#'
#' @param exp_mat A numeric matrix with features (e.g., probes, genes) as rows
#'   and samples as columns.
#'
#' @return A numeric matrix with the same dimensions as the input, but with rows
#'   reordered so that features with the highest standard deviation appear first.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Calculates the standard deviation for each row using \code{matrixStats::rowSds()}
#'   \item Determines the ordering of rows by decreasing standard deviation
#'   \item Reorders the matrix rows accordingly
#' }
#'
#' This approach is commonly used in genomics workflows to focus analyses on the
#' most informative (variable) features, as features with low variability contribute
#' little information to downstream analyses like PCA.
#'
#' @note Requires the \code{matrixStats} package for efficient standard deviation
#'   calculation on matrix rows.
OrderDataBySd <- function(exp_mat){
  # compute sds for each row
  sds <- matrixStats::rowSds(exp_mat)
  sdsSorted <- order(sds, decreasing = TRUE)
  
  # order by most variable probes on top
  exp_mat[sdsSorted ,]
}

#' Calculate Principal Component Analysis on DNA Methylation Data
#'
#' Performs principal component analysis (PCA) on DNA methylation beta values
#' after converting to M-values (logit transformation) and selecting the most
#' variable probes. Uses truncated SVD via irlba for computational efficiency
#' with large datasets.
#'
#' @param dnam_data A SummarizedExperiment object containing DNA methylation data
#'   with beta values in the "dnam" assay slot.
#' @param pc_count Integer. Number of principal components to compute. Default is 20.
#'
#' @return A prcomp_irlba object containing PCA results with the following elements:
#'   \describe{
#'     \item{x}{Matrix of principal component scores (samples x PCs)}
#'     \item{sdev}{Standard deviations of the principal components}
#'     \item{rotation}{Matrix of variable loadings (features x PCs)}
#'     \item{center}{Centers used for scaling}
#'     \item{scale}{Scales used for normalization}
#'     \item{totalvar}{Total variance in the data}
#'   }
#'
#' @details
#' The PCA workflow:
#' \enumerate{
#'   \item Extracts beta values from the SummarizedExperiment object
#'   \item Converts beta values to M-values using \code{minfi::logit2()} for
#'     better normal distribution properties
#'   \item Orders probes by standard deviation using \code{OrderDataBySd()}
#'   \item Selects the top 50,000 most variable probes
#'   \item Computes PCA on the transposed matrix (samples as rows) with centering
#'     and scaling
#'   \item Uses \code{irlba::prcomp_irlba()} for efficient computation of truncated SVD
#' }
#'
#' M-value transformation is preferred over beta values for PCA because:
#' \itemize{
#'   \item M-values have better statistical properties (more homoscedastic)
#'   \item Beta values are bounded [0,1], creating artificial constraints
#'   \item M-values better approximate normal distribution assumptions
#' }
#'
#' Using only the 50,000 most variable probes:
#' \itemize{
#'   \item Reduces computational burden
#'   \item Focuses on informative features
#'   \item Removes low-variance probes that contribute mainly noise
#' }
#'
#' @note
#' Requires:
#' \itemize{
#'   \item \code{irlba} package for efficient truncated SVD
#'   \item \code{minfi} package for logit2 transformation
#'   \item \code{matrixStats} package (via \code{OrderDataBySd()})
#' }
calculate_dnam_pca <- function(dnam_data, pc_count = 20) {
  beta_data <- SummarizedExperiment::assays(dnam_data)$dnam
  M_mat <- minfi::logit2(beta_data)
  M_mat <- OrderDataBySd(M_mat)
  pc_vals <- irlba::prcomp_irlba(
    t(M_mat[1:50000,]),
    n = pc_count,
    center = TRUE,
    scale = TRUE
  )
  
  pc_vals
}

#' Summarize DNA Methylation PCA Results and Identify Outliers
#'
#' Processes PCA results to calculate variance explained, extract PC scores,
#' compute z-scores, and identify outlier samples based on distance from the
#' mean in PC1 and PC2 space.
#'
#' @param pc_vals A prcomp_irlba object returned by \code{calculate_dnam_pca()}.
#' @param outlier_cutoff Numeric. Z-score threshold for identifying outliers.
#'   Samples with absolute z-score > outlier_cutoff on PC1 or PC2 are flagged.
#'   Common values are 3 (3 standard deviations) or 4 (more conservative).
#' @param pc_count Integer. Number of principal components to include in the
#'   output data frame. Default is 5.
#' @param pc_label Character. Prefix for principal component column names.
#'   Default is "PC_".
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{pc_df}{Data frame with PC scores (PC_1 to PC_N), z-scores for PC1 and PC2,
#'       outlier flags for each PC, and an overall is_outlier flag}
#'     \item{percentVar}{Numeric vector of proportion of variance explained by each PC}
#'     \item{pca_params}{Data frame with summary statistics (mean, sd, variance) for PC1 and PC2}
#'   }
#'
#' @details
#' The function performs several calculations:
#'
#' \strong{Variance Explained:}
#' Calculates the proportion of total variance explained by each PC using the
#' formula: \code{sdev^2 / totalvar}
#'
#' \strong{PCA Parameters:}
#' Computes summary statistics for PC1 and PC2:
#' \itemize{
#'   \item Mean scores
#'   \item Standard deviations
#'   \item Proportion of variance explained
#' }
#'
#' \strong{Outlier Detection:}
#' \enumerate{
#'   \item Calculates z-scores for each sample on PC1 and PC2
#'   \item Flags samples where \code{abs(z-score) > outlier_cutoff}
#'   \item Creates overall outlier flag (TRUE if outlier on either PC1 or PC2)
#' }
#'
#' Z-score calculation: \code{(score - mean) / sd}
#'
#' Outliers represent samples that are unusually distant from the population
#' center in PC space, potentially indicating:
#' \itemize{
#'   \item Technical issues (sample swaps, processing errors)
#'   \item Biological outliers (unusual methylation patterns)
#'   \item Batch effects
#'   \item Contamination
#' }
summarise_dnam_pca <- function(
    pc_vals, outlier_cutoff, pc_count = 5, pc_label = "PC_"
) {
  percentVar <- pc_vals$sdev^2 / pc_vals$totalvar
  
  pca_params <- data.frame(
    statistic = c(
      "pc1_mean", "pc2_mean", "pc1_sd", "pc2_sd", "pc1_var", "pc2_var"),
    value = c(
      mean(pc_vals$x[,1]),
      mean(pc_vals$x[,2]),
      stats::sd(pc_vals$x[,1]),
      stats::sd(pc_vals$x[,2]),
      percentVar[[1]],
      percentVar[[2]]
    )
  ) %>%
    tibble::column_to_rownames(var = "statistic")
  
  pc_df <- data.frame(pc_vals$x[,1:pc_count])
  colnames(pc_df) <- paste0(pc_label, 1:pc_count)
  pc1 <- as.numeric(pc_df[,paste0(pc_label, "1")])
  mean1 <- pca_params["pc1_mean", "value"]
  sd1 <- pca_params["pc1_sd", "value"]
  
  pc2 <- as.numeric(pc_df[,paste0(pc_label, "2")])
  mean2 <- pca_params["pc2_mean", "value"]
  sd2 <- pca_params["pc2_sd", "value"]
  
  pc_df$pc1_zscore <- (pc1 - mean1) / sd1
  pc_df$pc2_zscore <- (pc2 - mean2) / sd2
  pc_df <- pc_df %>%
    dplyr::mutate(
      pc1_outlier = abs(.data$pc1_zscore) > outlier_cutoff,
      pc2_outlier = abs(.data$pc2_zscore) > outlier_cutoff
    ) %>%
    dplyr::mutate(
      is_outlier = .data$pc1_outlier | .data$pc2_outlier
    )
  
  list(pc_df = pc_df, percentVar = percentVar, pca_params = pca_params)
}

#' Complete DNA Methylation PCA Workflow with Outlier Detection
#'
#' Performs principal component analysis on DNA methylation data, identifies
#' outlier samples, and integrates PC scores and outlier flags into the
#' SummarizedExperiment colData. This is the main function for PCA-based
#' quality control in DNA methylation studies.
#'
#' @param dnam_data A SummarizedExperiment object containing DNA methylation data
#'   with beta values in the "dnam" assay slot.
#' @param outlier_cutoff Numeric. Z-score threshold for identifying outliers.
#'   Samples with absolute z-score > outlier_cutoff on PC1 or PC2 are flagged
#'   as outliers. Common values: 3 (standard) or 4 (conservative).
#' @param pc_count Integer. Number of principal components to compute and include
#'   in the output. Default is 5.
#' @param pc_label Character. Prefix for principal component column names in colData.
#'   Default is "PC_".
#' @param outlier_label Character. Column name for the outlier flag in colData.
#'   Default is "is_outlier".
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{dnam_data}{Updated SummarizedExperiment with PC scores and outlier
#'       flags added to colData}
#'     \item{percentVar}{Numeric vector of proportion of variance explained by each PC}
#'     \item{pca_params}{Data frame with summary statistics for PC1 and PC2}
#'   }
#'
#' @details
#' This function combines the complete PCA workflow:
#' \enumerate{
#'   \item Computes PCA using \code{calculate_dnam_pca()}
#'   \item Summarizes results and identifies outliers using \code{summarise_dnam_pca()}
#'   \item Extracts existing phenotype data from colData
#'   \item Appends PC scores (PC_1 through PC_N) to phenotype data
#'   \item Adds outlier flag to phenotype data
#'   \item Reconstructs SummarizedExperiment with updated colData
#'   \item Returns both the updated object and PCA statistics
#' }
#'
#' The updated SummarizedExperiment colData will contain:
#' \itemize{
#'   \item Original phenotype columns
#'   \item New columns: PC_1, PC_2, ..., PC_N (or custom prefix)
#'   \item New column: is_outlier (or custom name) - logical flag
#' }
#'
#' This integrated approach ensures that PCA results are properly linked to
#' sample metadata and can be easily used in downstream analyses or plotting.
#'
#' @note
#' The function preserves all original assays and rowData, only updating colData.
get_dnam_pca <- function(
    dnam_data,
    outlier_cutoff,
    pc_count = 5,
    pc_label = "PC_",
    outlier_label = "is_outlier"
) {
  pc_vals <- calculate_dnam_pca(dnam_data, pc_count = pc_count)
  pca_res <- summarise_dnam_pca(
    pc_vals, outlier_cutoff, pc_count = pc_count, pc_label
  )
  pc_df <- pca_res$pc_df
  percentVar <- pca_res$percentVar
  pca_params <- pca_res$pca_params
  
  pheno_df <- dnam_data %>%
    SummarizedExperiment::colData() %>%
    as.data.frame()
  
  pheno_df <- cbind(pheno_df, pc_df[,paste0(pc_label, 1:pc_count)])
  pheno_df[,outlier_label] <- pc_df$is_outlier
  
  dnam_data <- SummarizedExperiment::SummarizedExperiment(
    assays = SummarizedExperiment::assays(dnam_data),
    rowData = SummarizedExperiment::rowData(dnam_data),
    colData = pheno_df
  )
  
  list(dnam_data = dnam_data, percentVar = percentVar, pca_params = pca_params)
}

#' Create PCA Scatter Plot with Outlier Highlighting
#'
#' Generates a ggplot2 scatter plot of PC1 vs PC2 with samples colored by a
#' specified variable, outliers labeled, and outlier boundaries marked with
#' dashed lines. This is useful for visualizing sample clustering and identifying
#' problematic samples.
#'
#' @param pca_df Data frame containing PCA results and sample metadata. Must
#'   include columns for sample IDs, PC scores, the plotting variable, and
#'   outlier flags.
#' @param dataset Character. Title for the plot (typically dataset name).
#' @param plot_var Character. Name of the column in \code{pca_df} to use for
#'   coloring points (e.g., "sex", "Batch", "tissue_type").
#' @param percentVar Numeric vector of variance proportions for each PC, used
#'   for axis labels.
#' @param pca_params Data frame with PCA summary statistics (means and SDs for
#'   PC1 and PC2), used to draw outlier boundaries.
#' @param outlier_cutoff Numeric. Z-score threshold used for outlier detection,
#'   used to position the outlier boundary lines.
#' @param pc_label Character. Prefix used for PC column names in \code{pca_df}.
#'   Default is "PC_".
#' @param outlier_label Character. Column name for outlier flag in \code{pca_df}.
#'   Default is "is_outlier".
#'
#' @return A ggplot2 object that can be displayed, saved, or further customized.
#'
#' @details
#' The plot includes several elements:
#'
#' \strong{Main Features:}
#' \itemize{
#'   \item Scatter plot of PC1 vs PC2
#'   \item Points colored by the specified variable
#'   \item Axis labels showing variance explained by each PC
#'   \item Title showing dataset name
#' }
#'
#' \strong{Outlier Visualization:}
#' \itemize{
#'   \item Dashed lines marking outlier boundaries (mean ± cutoff × SD)
#'   \item Text labels for outlier samples (using ggrepel to avoid overlaps)
#'   \item Four boundary lines: upper/lower for PC1 and PC2
#' }
#'
#' \strong{Column Requirements:}
#' The function validates that \code{pca_df} contains all required columns:
#' \itemize{
#'   \item "sample" - sample identifiers
#'   \item PC columns (e.g., "PC_1", "PC_2")
#'   \item \code{plot_var} - the variable for coloring
#'   \item \code{outlier_label} - outlier flags
#' }
#'
#' The boundary lines are calculated as:
#' \itemize{
#'   \item Upper boundary: mean + (cutoff × SD)
#'   \item Lower boundary: mean - (cutoff × SD)
#' }
plot_pca <- function(
    pca_df,
    dataset,
    plot_var,
    percentVar,
    pca_params,
    outlier_cutoff,
    pc_label = "PC_",
    outlier_label = "is_outlier"
) {
  columns <- c(
    "sample",
    paste0(pc_label, "1"),
    paste0(pc_label, "2"),
    plot_var,
    outlier_label
  )
  
  missing_col <- FALSE
  for (column in columns) {
    if (!column %in% colnames(pca_df)) {
      message("There is no '", column, "' column in the dataframe.")
      missing_col <- TRUE
    }
  }
  if (missing_col) {
    stop("At least 1 column is missing.")
  }
  
  column_names <- c("sample", "PC1", "PC2", plot_var, "is_outlier")
  
  plot_df <- pca_df[,columns]
  colnames(plot_df) <- column_names
  
  plot_df <- plot_df %>%
    dplyr::mutate(label_name = ifelse(.data$is_outlier, .data$sample, ""))
  subset_df <- plot_df %>%
    dplyr::filter(.data$is_outlier)
  
  p <- ggplot2::ggplot(
    data = plot_df, mapping = ggplot2::aes(
      x = PC1, y = PC2, color = !!sym(plot_var))) +
    ggplot2::geom_point(size = 1) +
    ggplot2::theme_bw() +
    ggplot2::xlab(paste0(
      "PC1: ", round(percentVar[[1]] * 100), "% variance")) +
    ggplot2::ylab(paste0(
      "PC2: ", round(percentVar[[2]] * 100), "% variance")) +
    ggplot2::geom_hline(
      yintercept = (pca_params["pc2_mean", "value"] +
                      outlier_cutoff * pca_params["pc2_sd", "value"]),
      linetype = "dashed") +
    ggplot2::geom_hline(
      yintercept = (pca_params["pc2_mean", "value"] -
                      outlier_cutoff * pca_params["pc2_sd", "value"]),
      linetype = "dashed") +
    ggplot2::geom_vline(
      xintercept = (pca_params["pc1_mean", "value"] +
                      outlier_cutoff * pca_params["pc1_sd", "value"]),
      linetype = "dashed") +
    ggplot2::geom_vline(
      xintercept = (pca_params["pc1_mean", "value"] -
                      outlier_cutoff * pca_params["pc1_sd", "value"]),
      linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = subset_df,
      ggplot2::aes(label = label_name),
      show.legend = FALSE, max.overlaps = 1000
    ) +
    ggplot2::ggtitle(dataset)
  
  p
}

#' Process Dataset with PCA and Optional Outlier Filtering
#'
#' Complete workflow function that performs PCA on DNA methylation data, identifies
#' outliers, and optionally filters outliers from the dataset. This is the main
#' entry point for PCA-based quality control with automatic sample removal.
#'
#' @param dnam_data A SummarizedExperiment object containing DNA methylation data
#'   with beta values in the "dnam" assay slot.
#' @param outlier_cutoff Numeric. Z-score threshold for identifying outliers.
#'   Default is 3 (3 standard deviations from the mean).
#' @param pc_count Integer. Number of principal components to compute. Default is 5.
#' @param pc_label Character. Prefix for PC column names in colData. Default is "PC_".
#' @param outlier_label Character. Column name for outlier flag in colData.
#'   Default is "is_outlier".
#' @param filter_outlier Logical. If TRUE, removes outlier samples from the returned
#'   dataset. If FALSE, keeps all samples but flags outliers. Default is TRUE.
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{dnam_data}{SummarizedExperiment object, filtered to remove outliers if
#'       \code{filter_outlier = TRUE}, otherwise containing all samples with outlier flags}
#'     \item{pca_df}{Data frame with complete sample metadata including PC scores and
#'       outlier flags for ALL samples (before filtering)}
#'     \item{percentVar}{Numeric vector of proportion of variance explained by each PC}
#'     \item{pca_params}{Data frame with summary statistics for PC1 and PC2}
#'   }
#'
#' @details
#' This function provides a complete PCA quality control workflow:
#'
#' \strong{Workflow Steps:}
#' \enumerate{
#'   \item Performs PCA and outlier detection using \code{get_dnam_pca()}
#'   \item Extracts updated metadata with PC scores and outlier flags
#'   \item If \code{filter_outlier = TRUE}: Subsets data to remove outlier samples
#'   \item If \code{filter_outlier = FALSE}: Keeps all samples with outlier flags
#'   \item Returns both filtered data and complete metadata
#' }
#'
#' \strong{Key Features:}
#' \itemize{
#'   \item \code{pca_df} always contains ALL samples (before filtering), allowing
#'     you to review which samples were removed
#'   \item \code{dnam_data} contains only non-outliers if \code{filter_outlier = TRUE}
#'   \item PC scores and outlier flags are preserved in colData
#'   \item Flexible outlier handling based on your QC strategy
#' }
#'
#' \strong{Typical Use Cases:}
#' \itemize{
#'   \item \code{filter_outlier = TRUE}: Automatic QC for production pipelines
#'   \item \code{filter_outlier = FALSE}: Exploratory analysis where you want to
#'     visualize outliers before deciding on removal
#' }
process_dataset <- function(
    dnam_data,
    outlier_cutoff = 3,
    pc_count = 5,
    pc_label = "PC_",
    outlier_label = "is_outlier",
    filter_outlier = TRUE
) {
  pca_res <- get_dnam_pca(
    dnam_data,
    outlier_cutoff,
    pc_count = pc_count, 
    pc_label = pc_label,
    outlier_label = outlier_label
  )
  dnam_data <- pca_res$dnam_data
  percentVar <- pca_res$percentVar
  pca_params <- pca_res$pca_params
  
  pca_df <- dnam_data %>%
    SummarizedExperiment::colData() %>%
    as.data.frame()
  
  if (filter_outlier) {
    dnam_data <- dnam_data[,pca_df[,outlier_label]]
  }
  
  list(
    dnam_data = dnam_data,
    pca_df = pca_df,
    percentVar = percentVar,
    pca_params = pca_params
  )
}