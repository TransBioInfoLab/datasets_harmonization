#' Perform Linear Regression for Single CpG Against Phenotype Parameter
#'
#' Fits a linear model to assess the association between a single CpG's
#' methylation values and a continuous phenotype parameter. Returns the
#' regression coefficient, t-value, and p-value for the CpG predictor.
#'
#' @param pheno_df Data frame containing phenotype data with samples as rows
#'   and phenotype variables as columns.
#' @param beta_value Numeric vector of DNA methylation beta values for a single
#'   CpG probe across samples. Length must match the number of rows in \code{pheno_df}.
#' @param parameter Character string specifying the column name in \code{pheno_df}
#'   to use as the response variable in the linear model.
#'
#' @return A one-row data frame with three columns:
#'   \describe{
#'     \item{estimate}{Regression coefficient (slope) for the CpG predictor}
#'     \item{t_value}{T-statistic for testing if the coefficient differs from zero}
#'     \item{p_value}{P-value from the t-test}
#'   }
#'
#' @details
#' The function fits a linear model of the form:
#' \code{parameter ~ cpg}
#'
#' Where:
#' \itemize{
#'   \item \code{parameter} is the phenotype response variable (e.g., age, BMI)
#'   \item \code{cpg} is the methylation beta value predictor
#' }
#'
#' This tests whether methylation at the CpG site is associated with the
#' phenotype parameter. A significant p-value indicates that methylation
#' levels correlate with the phenotype.
#'
#' The function is designed to be called repeatedly (e.g., via \code{apply()})
#' across many CpG probes to identify methylation sites associated with a
#' phenotype of interest.
run_lm_step <- function(pheno_df, beta_value, parameter) {
  formula <- stats::as.formula(paste0(parameter, " ~ cpg"))
  data <- data.frame(cpg = beta_value, pheno_df)
  lm_mod <- stats::lm(formula, data = data)
  coef_tbl <- summary(lm_mod)$coefficients %>%
    as.data.frame() %>%
    janitor::clean_names() %>%
    dplyr::select("estimate", "t_value", p_value = "pr_t")
  
  coef_tbl["cpg",]
}

#' Perform Logistic Regression for Single CpG Against Binary Phenotype
#'
#' Fits a logistic regression model to assess the association between a single
#' CpG's methylation values and a binary phenotype parameter. Returns the
#' regression coefficient, z-value, and p-value for the CpG predictor.
#'
#' @param pheno_df Data frame containing phenotype data with samples as rows
#'   and phenotype variables as columns.
#' @param beta_value Numeric vector of DNA methylation beta values for a single
#'   CpG probe across samples. Length must match the number of rows in \code{pheno_df}.
#' @param parameter Character string specifying the column name in \code{pheno_df}
#'   to use as the response variable. Should be a binary/dichotomous variable
#'   (e.g., case/control, diseased/healthy, 0/1).
#'
#' @return A one-row data frame with three columns:
#'   \describe{
#'     \item{estimate}{Regression coefficient (log odds ratio) for the CpG predictor}
#'     \item{z_value}{Z-statistic for testing if the coefficient differs from zero}
#'     \item{p_value}{P-value from the Wald test}
#'   }
#'
#' @details
#' The function fits a generalized linear model with binomial family:
#' \code{parameter ~ cpg}
#'
#' Where:
#' \itemize{
#'   \item \code{parameter} is the binary phenotype response variable
#'   \item \code{cpg} is the methylation beta value predictor
#' }
#'
#' The model uses a logit link function, testing whether methylation at the
#' CpG site is associated with the binary outcome. The coefficient estimate
#' represents the log odds ratio - the change in log odds of the outcome for
#' a unit increase in methylation.
#'
#' This approach is appropriate for:
#' \itemize{
#'   \item Case-control studies (disease vs. healthy)
#'   \item Binary traits (smoker vs. non-smoker)
#'   \item Dichotomized outcomes
#' }
#'
#' The function is designed to be called repeatedly across many CpG probes
#' to identify methylation sites associated with a binary phenotype.
run_glm_step <- function(pheno_df, beta_value, parameter) {
  formula <- stats::as.formula(paste0(parameter, " ~ cpg"))
  data <- data.frame(cpg = beta_value, pheno_df)
  glm_mod <- stats::glm(formula, data = data, family = "binomial")
  coef_tbl <- summary(glm_mod)$coefficients %>%
    as.data.frame() %>%
    janitor::clean_names() %>%
    dplyr::select("estimate", "z_value", p_value = "pr_z")
  
  coef_tbl["cpg",]
}

#' Run Statistical Tests Across All CpG Probes for Phenotype Association
#'
#' Performs genome-wide association testing between DNA methylation beta values
#' and a phenotype parameter across all CpG probes. Applies either linear
#' regression (continuous phenotype) or logistic regression (binary phenotype)
#' in parallel, then adjusts p-values for multiple testing using FDR correction.
#'
#' @param pheno_df Data frame containing phenotype data with samples as rows
#'   and phenotype variables as columns. Sample order must match columns in \code{beta}.
#' @param beta Numeric matrix of DNA methylation beta values with CpG probes as
#'   rows and samples as columns.
#' @param parameter Character string specifying the column name in \code{pheno_df}
#'   to test for association with methylation.
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'   Default is 16.
#' @param stat_test Character string specifying the statistical test: "lm" for
#'   linear regression (continuous phenotype) or "glm" for logistic regression
#'   (binary phenotype). Default is "lm".
#'
#' @return A data frame with one row per CpG probe containing:
#'   \describe{
#'     \item{probe}{CpG probe identifier}
#'     \item{estimate}{Regression coefficient}
#'     \item{t_value or z_value}{Test statistic (t for lm, z for glm)}
#'     \item{p_value}{Unadjusted p-value}
#'     \item{fdr}{FDR-adjusted p-value (Benjamini-Hochberg correction)}
#'   }
#'   Rows are ordered as in the input beta matrix.
#'
#' @details
#' The function performs an epigenome-wide association study (EWAS):
#'
#' \strong{Workflow:}
#' \enumerate{
#'   \item Validates the \code{stat_test} argument
#'   \item Selects appropriate test function (\code{run_lm_step} or \code{run_glm_step})
#'   \item Initializes parallel processing with specified cores
#'   \item Applies the test function to each row (CpG) of the beta matrix in parallel
#'   \item Combines results into a single data frame
#'   \item Stops parallel processing
#'   \item Adds FDR-adjusted p-values using Benjamini-Hochberg method
#' }
#'
#' \strong{Test Selection:}
#' \itemize{
#'   \item \strong{Linear regression (lm)}: For continuous phenotypes (age, BMI, etc.)
#'   \item \strong{Logistic regression (glm)}: For binary phenotypes (case/control, disease status)
#' }
#'
#' \strong{Multiple Testing Correction:}
#' FDR adjustment controls the expected proportion of false positives among
#' discoveries. An FDR < 0.05 means that among probes called significant,
#' we expect < 5% to be false positives.
run_stats_test <- function(
    pheno_df, beta, parameter, cores = 16, stat_test = c("lm", "glm")
) {
  stat_test <- match.arg(stat_test) 
  if (stat_test == "lm") {
    step_fun <- run_lm_step
  } else {
    step_fun <- run_glm_step
  }
  parallel_res <- start_parallel(TRUE, cores)
  do_parallel <- parallel_res$parallel
  cluster <- parallel_res$cluster
  
  result_df <- plyr::adply(
    beta,
    .margins = 1,
    .fun = function(beta_value) {
      step_fun(pheno_df, beta_value, parameter)
    }, .parallel = do_parallel
  )
  
  stop_parallel(do_parallel, cluster)
  
  result_df <- result_df %>%
    dplyr::rename(probe = "X1") %>%
    dplyr::mutate(
      fdr = stats::p.adjust(.data$p_value, method = "fdr")
    )
  
  result_df
}

#' Impute Missing Phenotype Data Using Methylation-Associated CpGs
#'
#' Performs informed imputation of missing phenotype values by first identifying
#' CpG probes significantly associated with the phenotype, then using those
#' probes' methylation values to guide multiple imputation with MICE.
#'
#' @param pheno_df Data frame containing phenotype data with samples as rows
#'   and phenotype variables as columns. The column specified by \code{parameter}
#'   may contain missing values (NA) to be imputed.
#' @param beta Numeric matrix of DNA methylation beta values with CpG probes as
#'   rows and samples as columns. Sample order must match rows in \code{pheno_df}.
#' @param parameter Character string specifying the column name in \code{pheno_df}
#'   containing the phenotype to impute.
#' @param fdr_thresh Numeric. FDR threshold for selecting significantly associated
#'   CpG probes. Only probes with FDR < \code{fdr_thresh} are used for imputation.
#'   Default is 0.01 (1% FDR).
#' @param stat_test Character string specifying the statistical test: "lm" for
#'   linear regression (continuous phenotype) or "glm" for logistic regression
#'   (binary phenotype). Default is "lm".
#'
#' @return A table showing the frequency of each value in the imputed phenotype
#'   variable. This allows verification that imputed values are reasonable.
#'
#' @details
#' This function implements an informed imputation strategy:
#'
#' \strong{Rationale:}
#' Standard imputation methods may not leverage the rich information in DNA
#' methylation data. By first identifying CpG probes associated with the
#' phenotype, we can use those methylation patterns to inform imputation,
#' potentially improving accuracy.
#'
#' \strong{Workflow:}
#' \enumerate{
#'   \item Validates the \code{stat_test} argument
#'   \item Performs EWAS using \code{run_stats_test()} to find phenotype-associated CpGs
#'   \item Filters to probes with FDR < \code{fdr_thresh}
#'   \item Extracts methylation data for significant probes
#'   \item Combines methylation data with the phenotype variable
#'   \item Performs multiple imputation using \code{mice::mice()} with:
#'     \itemize{
#'       \item m = 5 imputations
#'       \item seed = 42 for reproducibility
#'     }
#'   \item Returns the first imputed dataset (of 5 generated)
#'   \item Displays frequency table of imputed values
#' }
#'
#' \strong{Advantages:}
#' \itemize{
#'   \item Uses biologically relevant methylation patterns
#'   \item Focuses on probes with demonstrated association
#'   \item Multiple imputation accounts for uncertainty
#' }
#'
#' \strong{Limitations:}
#' \itemize{
#'   \item Requires sufficient non-missing data to identify associations
#'   \item May not work well if very few probes are significantly associated
#'   \item Assumes associations hold for samples with missing data
#' }
impute_parameter <- function(
    pheno_df, beta, parameter, fdr_thresh = 0.01, stat_test = c("lm", "glm")
) {
  stat_test <- match.arg(stat_test)
  result_df <- run_stats_test(
    pheno_df,
    beta,
    parameter, 
    cores = 16,
    stat_test = stat_test
  )
  
  probes <- result_df %>%
    dplyr::filter(.data$fdr < fdr_thresh) %>%
    dplyr::pull("probe")
  
  if (length(probes) < 50) {
    probes <- result_df %>%
      dplyr::arrange(.data$p_value) %>%
      head(n = 50) %>%
      dplyr::pull("probe")
    
    p_thresh <- result_df %>%
      dplyr::arrange(.data$p_value) %>%
      head(n = 50) %>%
      dplyr::pull("p_value") %>%
      max()
    
    message(
      "There were less than n = 50 probes at a threshold of FDR <= ",
      fdr_thresh,
      " so the n = 50 most significant probes were used at a threshold",
      " of p-value <= ",
      round(p_thresh, digits = 6),
      "."
    )
  } else {
    message(
      "At a threshold of FDR <= ", fdr_thresh,
      " n = ", length(probes),
      " probes were used for imputation."
    )    
  }
  
  if (stat_test == "glm") {
    method <- "logreg"
  } else {
    method <- "norm"
  }
  impute_data <- cbind(data.frame(t(beta[probes,])), test_var = pheno_df[,parameter])
  method <- c(rep("", length(probes)), method)
  imp <- mice::mice(impute_data, m = 5, seed = 42, method = method)
  res <- mice::complete(imp, 1)
  
  res$test_var
}
