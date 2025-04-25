#' RunEval
#'
#' @param fit A trained model object.
#' @param Test_set A data frame or matrix containing the test data. Rows represent samples, and columns represent features.
#' @param Test_label A data frame containing the test metadata. This should include columns for cohort identifiers cohortVar and class labels classVar.
#' @param Train_set A data frame or matrix containing the training data. Rows represent samples, and columns represent features.
#' @param Train_label A data frame containing the training metadata. This should include the class labels classVar.
#' @param Train_name A string specifying the cohort name for the training data. If not provided, defaults to "Training".
#' @param cohortVar A string specifying the column name in Test_label and Train_label that identifies cohorts
#' @param classVar A string specifying the column name in Test_label and Train_label that identifies the target class labels.
#'
#' @return
#' @export
#'
#' @examples
RunEval <- function(fit,
                    Test_set = NULL,
                    Test_label = NULL,
                    Train_set = NULL,
                    Train_label = NULL,
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    classVar){

  if(!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  }

  if((!is.null(Train_set)) & (!is.null(Train_label))) {
    new_data <- rbind.data.frame(Train_set[, fit$subFeature],
                                 Test_set[, fit$subFeature])

    if(!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Training"
    }
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    Test_label <- rbind.data.frame(Train_label[,c(cohortVar, classVar)],
                                   Test_label[,c(cohortVar, classVar)])
    Test_label[,1] <- factor(Test_label[,1],
                             levels = c(unique(Train_label[,cohortVar]), setdiff(unique(Test_label[,cohortVar]),unique(Train_label[,cohortVar]))))
  } else {
    new_data <- Test_set[, fit$subFeature]
  }

  RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))

  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS)
  Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
  unlist(lapply(Predict.out, function(data){
    as.numeric(auc(suppressMessages(roc(data[[classVar]], data$RS))))
  }))
}




