#' Start Parallel Processing Backend
#'
#' Initializes a parallel processing backend using doParallel and parallel packages.
#' Handles platform-specific differences between Windows and Unix-like systems.
#' On Windows, creates an explicit cluster object; on Unix-like systems, uses
#' forking via implicit parallelization.
#'
#' @param parallel Logical. Should parallel processing be enabled? If FALSE,
#'   no parallel backend is initialized regardless of other parameters.
#' @param cores Integer. The number of CPU cores to use for parallel processing.
#'   This value should typically not exceed the number of available cores on the system.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{parallel}{Logical indicating whether parallel processing was successfully initialized}
#'     \item{cluster}{A cluster object (Windows only) or NULL (Unix-like systems or if parallel=FALSE).
#'       This object is needed for cleanup with \code{stop_parallel()}}
#'   }
#'
#' @details
#' The function checks for the availability of required packages (doParallel and parallel)
#' before attempting to initialize parallel processing. If these packages are not available,
#' parallel processing is disabled and the function returns with parallel=FALSE.
#'
#' Platform-specific behavior:
#' \itemize{
#'   \item \strong{Windows}: Creates an explicit PSOCK cluster using \code{parallel::makeCluster()}
#'     and registers it with doParallel. Returns the cluster object for later cleanup.
#'   \item \strong{Unix-like (Linux/macOS)}: Uses implicit forking by calling
#'     \code{doParallel::registerDoParallel()} with the number of cores. Returns NULL
#'     for the cluster object as no explicit cleanup is needed.
#' }
#'
#'
#' @export
start_parallel <- function(parallel, cores) {
  if (parallel &&
      requireNamespace("doParallel", quietly = TRUE) &&
      requireNamespace("parallel", quietly = TRUE)) {
    if (Sys.info()["sysname"] == "Windows"){
      cluster <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cluster)
    } else {
      doParallel::registerDoParallel(cores)
      cluster <- NULL
    }
  } else {
    parallel = FALSE
    cluster <- NULL
  }
  
  list(parallel = parallel, cluster = cluster)
}

#' Stop Parallel Processing Backend
#'
#' Properly shuts down the parallel processing backend initialized by
#' \code{start_parallel()}. Handles cleanup of both implicit and explicit
#' cluster objects depending on the platform.
#'
#' @param parallel Logical. Indicates whether parallel processing was enabled.
#'   This should typically come from the return value of \code{start_parallel()}.
#' @param cluster A cluster object or NULL. On Windows systems, this is the
#'   explicit cluster object created by \code{start_parallel()}. On Unix-like
#'   systems, this should be NULL. This parameter should come from the return
#'   value of \code{start_parallel()}.
#'
#' @return Logical. Always returns TRUE indicating the function completed.
#'
#' @details
#' This function performs cleanup operations to properly shut down parallel
#' processing backends. It first stops the implicit cluster registered with
#' doParallel (if any), then on Windows systems, it stops the explicit cluster
#' object if provided.
#'
#' The function checks for the availability of required packages (doParallel
#' and parallel) before attempting cleanup operations. It is safe to call this
#' function even if parallel processing was not initialized, as it will simply
#' return TRUE without performing any operations.
stop_parallel <- function(parallel, cluster) {
  if (parallel &&
      requireNamespace("doParallel", quietly = TRUE) &&
      requireNamespace("parallel", quietly = TRUE)) {
    doParallel::stopImplicitCluster()
    if (!is.null(cluster) && Sys.info()["sysname"] == "Windows") {
      parallel::stopCluster(cluster)
    }
  }
  
  TRUE
}