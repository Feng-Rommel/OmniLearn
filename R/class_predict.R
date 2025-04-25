#' CalPredictScore
#'
#' @param fit A trained model object. Supported model types include "lognet", "glm", "svm.formula", "train", "glmboost", "plsRglmmodel", "rfsrc", "gbm", "xgb.Booster", and "naiveBayes".
#' @param new_data A data frame or matrix containing new data for prediction. Rows represent samples, and columns represent features.
#'
#' @return A named numeric vector of predicted values.
#' @export
#'
#' @examples
CalPredictScore <- function(fit, new_data){
  new_data <- new_data[, fit$subFeature]
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"      = predict(fit, type = 'response', as.matrix(new_data)), # response
    "glm"         = predict(fit, type = 'response', as.data.frame(new_data)), # response
    "svm.formula" = predict(fit, as.data.frame(new_data), probability = T), #
    "train"       = predict(fit, new_data, type = "prob")[[2]],
    "glmboost"    = predict(fit, type = "response", as.data.frame(new_data)), # response
    "plsRglmmodel" = predict(fit, type = "response", as.data.frame(new_data)), # response
    "rfsrc"        = predict(fit, as.data.frame(new_data))$predicted[, "1"],
    "gbm"          = predict(fit, type = 'response', as.data.frame(new_data)), # response
    "xgb.Booster" = predict(fit, as.matrix(new_data)),
    "naiveBayes" = predict(object = fit, type = "raw", newdata = new_data)[, "1"]
    # "drf" = predict(fit, functional = "mean", as.matrix(new_data))$mean
  ))
  RS = as.numeric(as.vector(RS))
  names(RS) = rownames(new_data)
  return(RS)
}


#' PredictClass
#'
#' @param fit A trained model object. Supported model types include "lognet", "glm", "svm.formula", "train", "glmboost", "plsRglmmodel", "rfsrc", "gbm", "xgb.Booster", and "naiveBayes"
#' @param new_data A data frame or matrix containing new data for prediction. Rows represent samples, and columns represent features.
#'
#' @return A named character vector of predicted class labels.
#' @export
#'
#' @examples
PredictClass <- function(fit, new_data){
  new_data <- new_data[, fit$subFeature]
  label <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"      = predict(fit, type = 'class', as.matrix(new_data)),
    "glm"         = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                           yes = "1", no = "0"), # glm不返回预测的类，将概率>0.5的作为1类
    "svm.formula" = predict(fit, as.data.frame(new_data), decision.values = T), #
    "train"       = predict(fit, new_data, type = "raw"),
    "glmboost"    = predict(fit, type = "class", as.data.frame(new_data)),
    "plsRglmmodel" = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                            yes = "1", no = "0"), # plsRglm不允许使用因子变量作为因变量，因而predict即使type设为class也无法正常运作
    "rfsrc"        = predict(fit, as.data.frame(new_data))$class,
    "gbm"          = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                            yes = "1", no = "0"), # gbm未设置预测类别，设置大于0.5为1
    "xgb.Booster" = ifelse(test = predict(fit, as.matrix(new_data))>0.5,
                           yes = "1", no = "0"), # xgboost 未提供预测类别，设置大于0.5为1
    "naiveBayes" = predict(object = fit, type = "class", newdata = new_data)
    # "drf" = predict(fit, functional = "mean", as.matrix(new_data))$mean
  ))
  label = as.character(as.vector(label))
  names(label) = rownames(new_data)
  return(label)
}
