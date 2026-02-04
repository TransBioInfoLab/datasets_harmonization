harmanScores_adjusted <- function(x, max_pc = 500) {
  x <- as.matrix(x)
  
  x <- scale(t(x), center = TRUE, scale = FALSE)
  cen <- attr(x, "scaled:center")
  sc <- attr(x, "scaled:scale")
  if(any(sc == 0)) {
    stop("cannot rescale a constant/zero column to unit variance")
  }
  
  s <- irlba::irlba(x, nv = min(nrow(x), max_pc))
  rotation <- s$v
  dimnames(rotation) <- list(
    colnames(x),
    paste0("PC", seq_len(ncol(rotation)))
  )
  scores <- x %*% rotation
  
  sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  
  r <- list(
    sdev = sdev,
    rotation = rotation,
    center = if(is.null(cen)) FALSE else cen,
    scale = if(is.null(sc)) FALSE else sc,
    scores = scores
  )
  
  class(r) <- "scores"
  r
}

harman_adjusted <- function(
    datamatrix,
    expt,
    batch,
    limit=0.95,
    numrepeats = 100000L,
    randseed,
    forceRand=FALSE,
    printInfo=FALSE,
    max_pc = 500
) {
  
  ######  Coerce a data.frame to a matrix  ##### 
  if(is.data.frame(datamatrix)) {
    datamatrix <- as.matrix(datamatrix)
  }
  
  ######  Sanity checks  #####
  if(!methods::is(datamatrix, "matrix")) {
    stop(paste("Require 'datamatrix' as a matrix or data.frame, not class \'",
               class(datamatrix), "\'.", sep=""))
  }
  
  if(!typeof(datamatrix) %in% c("integer", "double")) {
    stop(paste("'datamatrix' is type \'", typeof(datamatrix), "\',
               needs to be of type \'integer\' or \'double\'.", sep=""))
  }
  
  if(!is.vector(expt) && !is.factor(expt)) {
    stop(paste("Require 'expt' to be a vector or factor, not class \'",
               class(expt), "\'.", sep=""))
  }
  if(!is.vector(batch) && !is.factor(batch)) {
    stop(paste("Require 'batch' to be a vector or factor, not class \'",
               class(batch), "\'.", sep=""))
  }
  if(is.vector(expt)) {
    if(sum(is.na(expt)) > 0 ||
       sum(is.nan(expt)) > 0 ||
       sum(is.null(expt)) > 0) {
      stop("Cannot have NA, NaN or NULL as 'expt' levels.")
    }
  }
  if(is.vector(batch)) {
    if(sum(is.na(batch)) > 0 ||
       sum(is.nan(batch)) > 0 ||
       sum(is.null(batch)) > 0) {
      stop("Cannot have NA, NaN or NULL as 'batch' levels.")
    }
  }
  if(length(expt) != length(batch)) stop("'expt' and 'batch' vectors not the
                                         same length.")
  if(!is.numeric(limit) || limit < 0 || limit > 1) {
    stop(paste("'limit' needs to be a number between 0 and 1, not \"", limit,
               "\".", sep=""))
  }
  if(!is.numeric(numrepeats)) stop("'numrepeats' needs to be numeric.")
  
  if (!missing(randseed)) {
    if(!is.numeric(randseed)) stop("'randseed' needs to be numeric.")
  } else {
    randseed <- stats::runif(1, 0, 1e9)
  }
  
  strict <- FALSE
  #  Sanity checking to see if the expt vector length is equal to the
  #  number of matrix columns
  if(length(expt) != ncol(datamatrix)) {
    msg <- "'expt' vector not equal to the number of datamatrix columns."
    if(strict == FALSE) {
      warning(msg)
    } else {
      stop(msg)
    }
  }
  
  # Coerce expt and batch to factors
  expt <- factor(expt)
  batch <- factor(batch)
  
  if(length(levels(expt)) < 2 || length(levels(batch)) < 2) {
    stop("Require more than one experimental factor and/or batch for experiment
         structure")
  }
  
  #####  PCA  #####  
  
  # Don't shift the original data into the .RunHarman function as we just need
  # the PCs to kick it off.
  if(printInfo == TRUE) cat('Performing PCA... ')
  pca <- harmanScores_adjusted(datamatrix, max_pc = max_pc)
  pc_data_scores <- pca$scores
  # Try and free up RAM, but keep the sample names first.
  sample_names <- dimnames(datamatrix)[[2]]
  rm(datamatrix)
  gc()
  if(printInfo == TRUE) cat('done.\n')
  
  #####  Call Rcpp wrapper function  #####
  
  # Form group construct by converting all expt and batch names to an integer
  group <- as.matrix(data.frame(expt=as.integer(expt), batch=as.integer(batch)))
  rownames(group) <- sample_names
  
  
  if(printInfo == TRUE) cat('Now calling the Rcpp layer.\n')
  res <- Harman:::.callHarman(
    pc_data_scores,
    group,
    limit,
    numrepeats,
    randseed,
    forceRand,
    printInfo
  )
  
  #####  Form S3 object  #####
  
  parameters <- list(limit=limit, numrepeats=numrepeats, randseed=randseed,
                     forceRand=forceRand)
  factors <- data.frame(expt=expt, batch=batch)
  rownames(factors) <- sample_names
  factors$expt.numeric <- group[, 'expt']
  factors$batch.numeric <- group[, 'batch']
  dim_names <- paste('PC',seq_len(length(res[["confidence_vector"]])), sep='')
  stats <- data.frame(dimension=dim_names,
                      confidence=res[["confidence_vector"]],
                      correction=res[["correction_vector"]])
  
  # Corrected scores are returned as a numeric array
  # Need to coerce them into a matrix
  corrected <- matrix(res[["corrected_scores"]],
                      nrow=nrow(pc_data_scores),
                      ncol=ncol(pc_data_scores),
                      dimnames=dimnames(pc_data_scores))
  
  results <- list(factors=factors,
                  parameters=parameters,
                  stats=stats,
                  center=pca$center,
                  rotation=pca$rotation,
                  original=pc_data_scores,
                  corrected=corrected)
  
  # Define the S3 results class
  results <- structure(results, class = "harmanresults")
  results
}

reconstructData_adjusted <- function(object, this='corrected')  {
  
  # The Matlab code from which this function is built:
  # reconstructedmicroarray = corrected_scorebatch_matrix*coeff'+\
  # ones(n,1)*means_initialprobesets;
  # corrected_scorebatch_matrix == object[[this]]
  # coeff == t(object$rotation), so coeff' == object$rotation
  # n == number of samples
  # ones(n,1) == matrix(1, n, 1)
  # means_initialprobesets == object$center
  
  if(!methods::is(object, "harmanresults")) {
    stop(paste("Require an instance of 'harmanresults', not class \'",
               class(object), "\'.", sep=""))
  }
  
  if(!(this %in% c('original', 'corrected'))) {
    stop("Require 'this' to be either the values 'original' or 'corrected'.")
  }
  
  # Add the extra column of zeros.
  scores <- object[[this]]
  n <- nrow(scores)
  ones <- matrix(1, n, 1)
  
  # corrected_scorebatch_matrix*coeff'
  centred_scores_matrix <- scores %*% t(object$rotation)
  # ones(n,1)*means_initialprobesets
  means_matrix <- ones %*% matrix(object$center, nrow=1)
  # Now add the means and transpose
  reconstructed <- t(centred_scores_matrix + means_matrix)
  reconstructed
}