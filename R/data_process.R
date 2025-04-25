#' Standardize Data by Cohorts
#'
#' @param data A data frame or matrix where rows represent samples and columns represent genes.
#' The dataset to be standardized.
#' @param cohort An optional vector indicating the cohort (group) assignment for each sample.
#' If not provided (NULL), all samples are assumed to belong to a single cohort.
#' @param centerFlags Logical value(s) or named logical vector. Determines whether to center
#' the gene values to have a mean of zero within each cohort. Defaults to NULL (no centering).
#' If a single TRUE/FALSE value is provided, it applies to all cohorts. If a named vector is
#' provided, the names must correspond to cohort names.
#' @param scaleFlags Logical value(s) or named logical vector. Determines whether to scale
#' the gene values to have a standard deviation of one within each cohort. Defaults to NULL
#' (no scaling). The usage is similar to centerFlags.
#'
#' @return A standardized data matrix with rows reordered to match the original sample names.
#' @export
#'
#' @examples
#' data <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' cohort <- rep(c("A", "B"), each = 5)
#' centerFlags <- c("A" = TRUE, "B" = FALSE)
#' scaleFlags <- TRUE
#' scaledData <- scaleData(data, cohort, centerFlags, scaleFlags)
scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)
  if (is.null(cohort)){
    data <- list(data); names(data) = "training"
  }else{
    data <- split(as.data.frame(data), cohort)
  }

  if (is.null(centerFlags)){
    centerFlags = F; message("No centerFlags found, set as FALSE")
  }
  if (length(centerFlags)==1){
    centerFlags = rep(centerFlags, length(data)); message("set centerFlags for all cohort as ", unique(centerFlags))
  }
  if (is.null(names(centerFlags))){
    names(centerFlags) <- names(data); message("match centerFlags with cohort by order\n")
  }

  if (is.null(scaleFlags)){
    scaleFlags = F; message("No scaleFlags found, set as FALSE")
  }
  if (length(scaleFlags)==1){
    scaleFlags = rep(scaleFlags, length(data)); message("set scaleFlags for all cohort as ", unique(scaleFlags))
  }
  if (is.null(names(scaleFlags))){
    names(scaleFlags) <- names(data); message("match scaleFlags with cohort by order\n")
  }

  centerFlags <- centerFlags[names(data)]; scaleFlags <- scaleFlags[names(data)]
  outdata <- mapply(standarize.fun, indata = data, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = F)
  # lapply(out.data, function(x) summary(apply(x, 2, var)))
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  return(outdata)
}
