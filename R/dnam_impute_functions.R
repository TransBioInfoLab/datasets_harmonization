#' Extract beta-value matrix and apply detection p-value masking
#' Converts a DNAm SummarizedExperiment into a sample-by-probe beta matrix
#' and replaces beta values with NA for probe–sample pairs that failed
#' detection p-value thresholds.
#' @param dnam_data A SummarizedExperiment containing DNAm beta values
#'   in the \code{assays(dnam_data)$dnam} slot.
#' @param detp_df A data frame with columns \code{sample} and \code{cpg}
#'   indicating probe–sample combinations failing detection p-value QC.
#' @return A numeric matrix of beta values with samples in rows and probes
#'   in columns, with failed measurements set to NA.
#' @importFrom SummarizedExperiment assays
#' @importFrom dplyr filter pull
get_beta_data <- function(dnam_data, detp_df) {
  beta_data <- SummarizedExperiment::assays(dnam_data)$dnam %>%
    t()
  
  samples <- row.names(beta_data)
  probes <- colnames(beta_data)
  
  detp_df <- detp_df %>%
    dplyr::filter(.data$sample %in% samples) %>%
    dplyr::filter(.data$cpg %in% probes)
  
  for (Sample in unique(detp_df$sample)) {
    cpgs <- detp_df %>%
      dplyr::filter(.data$sample == Sample) %>%
      dplyr::pull("cpg")
    
    beta_data[Sample, cpgs] <- NA
  }
  
  beta_data
}

#' Extract probe annotation for imputation
#' Retrieves minimal probe-level annotation required for methylation
#' imputation, including CpG probe IDs and chromosome assignments.
#' @param dnam_data A SummarizedExperiment containing probe annotation
#'   in its \code{rowData}.
#' @return A data frame with columns \code{cpg} (probe ID) and \code{chr}
#'   (chromosome).
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select
get_probe_df <- function(dnam_data) {
  probe_df <- dnam_data %>%
    SummarizedExperiment::rowData() %>%
    as.data.frame() %>%
    dplyr::select(cpg = "Name", chr = "seqnames")
  
  probe_df
}

#' Impute missing beta values using batch-aware imputation
#' Performs methylation beta-value imputation using \pkg{methyLImp2},
#' incorporating batch structure to constrain imputation within batches.
#' @param beta_data A numeric matrix of beta values with samples in rows
#'   and probes in columns.
#' @param probe_df A data frame containing probe annotations required
#'   by \pkg{methyLImp2}.
#' @param batches A vector indicating batch membership for each sample,
#'   in the same order as rows of \code{beta_data}.
#' @return A numeric matrix of imputed beta values with original row and
#'   column names restored.
#' @importFrom BiocParallel SnowParam
#' @importFrom methyLImp2 methyLImp2
impute_beta <- function(beta_data, probe_df, batches) {
  samples <- row.names(beta_data)
  probes <- colnames(beta_data)
  
  beta_data <- methyLImp2::methyLImp2(
    beta_data,
    type = "user",
    annotation = probe_df,
    BPPARAM = BiocParallel::SnowParam(workers = 22, exportglobals = FALSE),
    groups = batches,
    overwrite_res = FALSE
  )
  
  colnames(beta_data) <- probes
  row.names(beta_data) <- samples
  
  beta_data
}

#' Impute missing beta values without batch constraints
#' Performs methylation beta-value imputation using \pkg{methyLImp2}
#' without incorporating batch structure, using minibatch-based
#' stochastic imputation.
#' @param beta_data A numeric matrix of beta values with samples in rows
#'   and probes in columns.
#' @param probe_df A data frame containing probe annotations required
#'   by \pkg{methyLImp2}.
#' @return A numeric matrix of imputed beta values with original row and
#'   column names restored.
#' @importFrom BiocParallel SnowParam
#' @importFrom methyLImp2 methyLImp2
impute_beta_nobatch <- function(beta_data, probe_df) {
  samples <- row.names(beta_data)
  probes <- colnames(beta_data)
  
  beta_data <- methyLImp2::methyLImp2(
    beta_data,
    type = "user",
    annotation = probe_df,
    BPPARAM = BiocParallel::SnowParam(workers = 22, exportglobals = FALSE),
    groups = NULL,
    overwrite_res = FALSE,
    minibatch_frac = 0.2,
    minibatch_reps = 3
  )
  
  colnames(beta_data) <- probes
  row.names(beta_data) <- samples
  
  beta_data
}

#' Impute beta values for a subset of batches
#' Subsets a DNAm SummarizedExperiment to selected batches, applies
#' detection p-value masking, performs batch-aware beta-value imputation,
#' and reconstructs a SummarizedExperiment with imputed values.
#' @param dnam_data A SummarizedExperiment containing DNAm beta values
#'   and batch information.
#' @param detp_df A data frame of detection p-value failures.
#' @param batches A character vector of batch IDs to include.
#' @return A SummarizedExperiment containing imputed beta values for the
#'   specified batches.
#' @importFrom SummarizedExperiment assays rowData colData SummarizedExperiment
#' @importFrom dplyr pull
impute_batches <- function(dnam_data, detp_df, batches) {
  batch_ls <- dnam_data %>%
    SummarizedExperiment::colData() %>%
    as.data.frame() %>%
    dplyr::pull("Batch")
  
  dnam_data <- dnam_data[,batch_ls %in% batches]
  beta_data <- get_beta_data(dnam_data, detp_df)
  probe_df <- get_probe_df(dnam_data)
  batches <- dnam_data %>%
    SummarizedExperiment::colData() %>%
    as.data.frame() %>%
    dplyr::pull("Batch")
  beta_data <- impute_beta(beta_data, probe_df, batches)
  
  dnam_data <- SummarizedExperiment::SummarizedExperiment(
    assays = list(dnam = t(beta_data)),
    rowData = dnam_data %>%
      SummarizedExperiment::rowData() %>%
      as.data.frame(),
    colData = dnam_data %>%
      SummarizedExperiment::colData() %>%
      as.data.frame()
  )
  
  dnam_data
}


#' Summarize probe missingness by batch
#' Computes per-batch summaries of probe-level missingness based on
#' detection p-value failures.
#' @param detp_df A data frame containing columns \code{cpg} and
#'   \code{Batch} indicating detection p-value failures.
#' @param pheno_df A phenotype data frame containing batch information.
#' @return A data frame summarizing total, present, and missing probe
#'   counts per batch and CpG.
#' @importFrom dplyr group_by summarise ungroup left_join mutate select
get_batch_summaries <- function(detp_df, pheno_df) {
  missing_df <- detp_df %>%
    dplyr::select("cpg", "Batch") %>%
    dplyr::group_by(.data$Batch, .data$cpg) %>%
    dplyr::summarise(missing_count = n()) %>%
    dplyr::ungroup()
  
  count_df <- pheno_df %>%
    dplyr::group_by(.data$Batch) %>%
    dplyr::summarise(total_count = n()) %>%
    dplyr::ungroup()
  
  batch_df <- missing_df %>%
    dplyr::left_join(count_df, by = "Batch") %>%
    dplyr::mutate(
      present_count = .data$total_count - .data$missing_count,
      label = paste0(.data$Batch, "_", .data$cpg)
    ) %>%
    dplyr::select(
      "cpg", "Batch", "label", "total_count", "present_count", "missing_count"
    )
  
  batch_df
}

#' Identify probe–sample pairs targeted for imputation
#' Identifies CpG–sample combinations to be imputed based on low probe
#' presence within batches, using detection p-value failure patterns.
#' @param pheno_df A phenotype data frame containing \code{sample} and
#'   \code{Batch} columns.
#' @param detp_df A data frame of detection p-value failures with
#'   \code{sample} and \code{cpg} columns.
#' @return A data frame containing \code{cpg} and \code{sample} columns
#'   specifying imputation targets.
#' @importFrom dplyr left_join filter select mutate
get_na_targets <- function(pheno_df, detp_df) {
  detp_df <- detp_df %>%
    dplyr::left_join(
      pheno_df %>% dplyr::select("sample", "Batch"),
      by = "sample"
    ) %>%
    dplyr::select("cpg", "sample", "Batch") %>%
    dplyr::mutate(label = paste0(.data$Batch, "_", .data$cpg))
  
  print(dim(detp_df))
  
  batch_df <- get_batch_summaries(detp_df, pheno_df) %>%
    dplyr::filter(.data$present_count < 10)
  
  detp_df <- detp_df %>%
    dplyr::filter(.data$label %in% batch_df$label) %>%
    dplyr::select("cpg", "sample")
  
  print(dim(detp_df))
  
  detp_df
}

#' Impute missing beta values without batch structure
#'
#' Performs methylation beta-value imputation using \pkg{methyLImp2}
#' without incorporating batch information. Imputation is carried out
#' using stochastic minibatching to improve robustness when batch
#' labels are unavailable or intentionally ignored.
#'
#' @param beta_data A numeric matrix of beta values with samples in rows
#'   and CpG probes in columns. Missing values (NA) indicate positions
#'   to be imputed.
#' @param probe_df A data frame containing probe-level annotation required
#'   by \pkg{methyLImp2}, typically including CpG probe IDs and chromosome
#'   assignments.
#'
#' @return A numeric matrix of beta values with imputed entries, preserving
#'   the original sample (row) and probe (column) names.
#'
#' @details
#' This function uses minibatch-based imputation (\code{minibatch_frac = 0.2},
#' \code{minibatch_reps = 3}) to stabilize estimates in the absence of
#' batch constraints. Parallel computation is enabled via
#' \code{BiocParallel::SnowParam}.
#'
#' @importFrom methyLImp2 methyLImp2
#' @importFrom BiocParallel SnowParam
impute_beta_nobatch <- function(beta_data, probe_df) {
  samples <- row.names(beta_data)
  probes <- colnames(beta_data)
  
  beta_data <- methyLImp2::methyLImp2(
    beta_data,
    type = "user",
    annotation = probe_df,
    BPPARAM = BiocParallel::SnowParam(workers = 22, exportglobals = FALSE),
    groups = NULL,
    overwrite_res = FALSE,
    minibatch_frac = 0.2,
    minibatch_reps = 3
  )
  
  colnames(beta_data) <- probes
  row.names(beta_data) <- samples
  
  beta_data
}