#' Extract Batch and Experimental Variables from DNA Methylation Data
#'
#' Extracts batch identifiers and experimental variables from a SummarizedExperiment
#' object containing DNA methylation data. Combines multiple experimental variables
#' into a single factor by concatenating them with underscores.
#'
#' @param dnam_data A SummarizedExperiment object containing DNA methylation data
#'   with column metadata (colData) that includes batch and experimental variables.
#' @param batch_var Character string specifying the column name in colData that
#'   contains batch identifiers. Default is "Batch".
#' @param expt_vars Character vector specifying one or more column names in colData
#'   that contain experimental variables (e.g., biological covariates). Default is
#'   \code{c("sex", "DX")}. If multiple variables are provided, they are concatenated
#'   with underscores to create a single experimental factor.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{batch}{A factor vector of batch identifiers, one per sample}
#'     \item{expt}{A factor vector of experimental conditions, one per sample.
#'       If multiple experimental variables are provided, they are combined as
#'       "var1_var2_var3"}
#'   }
#'
#' @details
#' This function prepares batch and experimental variables for batch correction
#' algorithms that require these inputs as factors. The experimental variables
#' represent biological or technical covariates that should be preserved during
#' batch correction.
#'
#' When multiple experimental variables are provided:
#' \itemize{
#'   \item They are extracted from colData in the order specified
#'   \item Each variable is converted to character
#'   \item Variables are concatenated with underscore separators
#'   \item The final combined string is converted to a factor
#' }
#'
#' For example, if \code{expt_vars = c("sex", "DX")} and a sample has
#' sex="M" and DX="Control", the resulting experimental factor level will be "M_Control".
#'
get_batch_params <- function(
    dnam_data, batch_var = "Batch", expt_vars = c("sex", "DX")
) {
  pheno_df <- dnam_data %>%
    SummarizedExperiment::colData() %>%
    as.data.frame()
  
  batches <- pheno_df %>%
    dplyr::pull(as.name(batch_var)) %>%
    as.character()
  
  expt <- pheno_df %>%
    dplyr::pull(as.name(expt_vars[[1]])) %>%
    as.character()
  
  if (length(expt_vars) > 1) {
    for (expt_var in expt_vars[2:length(expt_vars)]) {
      expt_add <- pheno_df %>%
        dplyr::pull(as.name(expt_var)) %>%
        as.character()
      expt <- paste0(expt, "_", expt_add)
    }
  }
  
  expt <- as.factor(expt)
  batches <- as.factor(batches)
  
  list(batch = batches, expt = expt)
}

#' Perform Batch Correction on DNA Methylation Data Using Harman
#'
#' Applies Harman batch correction to DNA methylation beta values after converting
#' them to M-values. Handles both full and adjusted (reduced dimensionality) Harman
#' approaches depending on the number of samples and specified maximum principal
#' components.
#'
#' @param dnam_data A SummarizedExperiment object containing DNA methylation data
#'   with beta values in the "dnam" assay and relevant metadata in colData.
#' @param batch_var Character string specifying the column name in colData that
#'   contains batch identifiers. Default is "Batch".
#' @param expt_vars Character vector specifying one or more column names in colData
#'   that contain experimental variables to preserve during batch correction.
#'   Default is \code{c("sex", "DX")}.
#' @param randseed Integer. Random seed for reproducibility of the Harman algorithm.
#'   Default is 42.
#' @param max_pc Integer. Maximum number of principal components to use in the
#'   adjusted Harman approach. If \code{max_pc >= ncol(dnam_data)}, the full
#'   Harman method is used instead. Default is 1000.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{harman_data}{The Harman object containing batch correction results.
#'       This can be either a standard Harman object (if \code{use_full = TRUE})
#'       or an adjusted Harman object (if \code{use_full = FALSE})}
#'     \item{use_full}{Logical indicating whether the full Harman method was used
#'       (\code{TRUE}) or the adjusted method with reduced dimensionality (\code{FALSE}).
#'       This flag is needed for proper data reconstruction}
#'   }
#'
#' @details
#' The function performs batch correction in several steps:
#' \enumerate{
#'   \item Extracts batch and experimental variables using \code{get_batch_params()}
#'   \item Converts beta values to M-values using \code{lumi::beta2m()}
#'   \item Applies a small shift (1e-4) to beta values before conversion to avoid
#'     extreme M-values from beta values near 0 or 1
#'   \item Determines whether to use full or adjusted Harman based on \code{max_pc}
#'     relative to sample size
#'   \item Runs Harman with a variance limit of 0.95 (95% of variance explained)
#' }
#'
#' \strong{Full vs Adjusted Harman:}
#' \itemize{
#'   \item \strong{Full}: Uses \code{Harman::harman()} when \code{max_pc >= ncol(dnam_data)}.
#'     Appropriate for datasets with moderate sample sizes.
#'   \item \strong{Adjusted}: Uses \code{harman_adjusted()} when sample size exceeds
#'     \code{max_pc}. Reduces computational burden by limiting principal components
#'     to \code{max_pc}, suitable for large datasets.
#' }
batch_data <- function(
    dnam_data,
    batch_var = "Batch",
    expt_vars = c("sex", "DX"),
    randseed = 42,
    max_pc = 1000
) {
  batch_params <- get_batch_params(
    dnam_data,
    batch_var = batch_var,
    expt_vars = expt_vars
  )
  batch <- batch_params$batch
  expt <- batch_params$expt
  
  m_data <- SummarizedExperiment::assays(dnam_data)$dnam %>%
    Harman::shiftBetas(shiftBy = 1e-4) %>%
    lumi::beta2m()
  
  if (max_pc >= ncol(dnam_data)) {
    use_full <- TRUE
    harman_data <- Harman::harman(
      m_data,
      expt = expt,
      batch = batch,
      limit = 0.95,
      randseed = randseed
    )    
  } else {
    use_full <- FALSE
    harman_data <- harman_adjusted(
      m_data,
      expt = expt,
      batch = batch,
      limit = 0.95,
      randseed = randseed,
      max_pc = max_pc
    )
  }
  
  list(harman_data = harman_data, use_full = use_full)
}

#' Reconstruct Batch-Corrected DNA Methylation Data
#'
#' Reconstructs batch-corrected DNA methylation data from Harman results, converts
#' M-values back to beta values, and packages the corrected data into a
#' SummarizedExperiment object with original metadata preserved.
#'
#' @param dnam_data The original SummarizedExperiment object containing DNA methylation
#'   data. This is used to extract and preserve the original row and column metadata.
#' @param harman_data A Harman object containing batch correction results, as returned
#'   by \code{batch_data()}. Can be either a full or adjusted Harman object.
#' @param use_full Logical indicating whether the full Harman method was used
#'   (\code{TRUE}) or the adjusted method (\code{FALSE}). This determines which
#'   reconstruction function to use. Default is \code{TRUE}. This value should
#'   come from the \code{batch_data()} return value.
#'
#' @return A SummarizedExperiment object containing:
#'   \describe{
#'     \item{assays}{A "dnam" assay with batch-corrected beta values}
#'     \item{rowData}{Original probe annotations from input \code{dnam_data}}
#'     \item{colData}{Original sample metadata from input \code{dnam_data}}
#'   }
#'
#' @details
#' The reconstruction process involves:
#' \enumerate{
#'   \item Reconstruct corrected M-values from Harman results
#'     \itemize{
#'       \item If \code{use_full = TRUE}: Uses \code{Harman::reconstructData()}
#'       \item If \code{use_full = FALSE}: Uses \code{reconstructData_adjusted()}
#'     }
#'   \item Convert corrected M-values back to beta values using \code{lumi::m2beta()}
#'   \item Package into SummarizedExperiment with original metadata preserved
#' }
#'
#' The function ensures that all probe annotations (rowData) and sample metadata
#' (colData) from the original dataset are preserved in the corrected dataset,
#' maintaining compatibility with downstream analyses.
reconstruct_data <- function(dnam_data, harman_data, use_full = TRUE) {
  if (use_full) {
    correct_data <- Harman::reconstructData(harman_data)
  } else {
    correct_data <- reconstructData_adjusted(harman_data)
  }
  
  correct_data <- lumi::m2beta(correct_data)
  correct_data <- SummarizedExperiment::SummarizedExperiment(
    assays = list(dnam = correct_data),
    rowData = dnam_data %>%
      SummarizedExperiment::rowData(),
    colData = dnam_data %>%
      SummarizedExperiment::colData()
  )
  
  correct_data
}