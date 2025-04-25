#' RunML
#'
#' @param method A string specifying the method used to create the model. Must be a supported method.
#' @param Train_set A matrix representing the training data. Rows correspond to samples, and columns correspond to genes.
#' @param Train_label Rows correspond to samples, and columns represent phenotypes. The phenotype should be encoded as 0 or 1.
#' @param mode A string with two possible values: 'Model' (default) or 'Variable'. 'Model' trains a model, while 'Variable' filters features.
#' @param classVar A string specifying the column name in Train_label to use as the target phenotype.
#'
#' @return A trained model or filtered features, depending on the selected mode.
#' @export
#'
#' @examples
RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
  # for example: Enet [alpha=0.4]
  method = gsub(" ", "", method) # 去除参数中的空格，得到Enet [alpha=0.4]
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)  # get name of ML algorithm, e.g., Enet
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method) # get parameter of ML algorithm, e.g., alpha=0.4

  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL
  )
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_set), " Variables")

  args = list("Train_set" = Train_set,
              "Train_label" = Train_label,
              "mode" = mode,
              "classVar" = classVar)
  args = c(args, method_param)

  obj <- do.call(what = paste0("Run", method_name),
                 args = args)

  if(mode == "Variable"){
    message(length(obj), " Variables retained;\n")
  }else{message("\n")}
  return(obj)
}

RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  cv.fit = cv.glmnet(x = Train_set,
                     y = Train_label[[classVar]],
                     family = "binomial", alpha = alpha, nfolds = 5)
  fit = glmnet(x = Train_set,
               y = Train_label[[classVar]],
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLasso <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
}

RunRidge <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
}

RunStepglm <- function(Train_set, Train_label, mode, classVar, direction){
  fit <- step(glm(formula = Train_label[[classVar]] ~ .,
                  family = "binomial",
                  data = as.data.frame(Train_set)),
              direction = direction, trace = 0)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunSVM <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit = svm(formula = eval(parse(text = paste(classVar, "~."))),
            data= data, probability = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLDA <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit = train(eval(parse(text = paste(classVar, "~."))),
              data = data,
              method="lda",
              trControl = trainControl(method = "cv"))
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunglmBoost <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial())

  cvm <- cvrisk(fit, papply = lapply,
                folds = cv(model.weights(fit), type = "kfold"))
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial(),
                  control = boost_control(mstop = max(mstop(cvm), 40)))

  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunplsRglm <- function(Train_set, Train_label, mode, classVar){
  cv.plsRglm.res = cv.plsRglm(formula = Train_label[[classVar]] ~ .,
                              data = as.data.frame(Train_set),
                              nt=min(5, ncol(Train_set)), verbose = F,
                              K=41)
  fit <- plsRglm(Train_label[[classVar]],
                 as.data.frame(Train_set),
                 modele = "pls-glm-logistic",
                 verbose = F, sparse = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunRF <- function(Train_set, Train_label, mode, classVar){
  rf_nodesize = 5 # may modify
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  fit <- rfsrc(formula = formula(paste0(classVar, "~.")),
               data = cbind(Train_set, Train_label[classVar]),
               ntree = 1000, nodesize = rf_nodesize,
               importance = T,
               proximity = T,
               forest = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunGBM <- function(Train_set, Train_label, mode, classVar){
  fit <- gbm(formula = Train_label[[classVar]] ~ .,
             data = as.data.frame(Train_set),
             distribution = 'bernoulli',
             n.trees = 5000,
             interaction.depth = 3,
             n.minobsinnode = 3,
             bag.fraction = 0.8,  # 减小袋装抽样比例
             shrinkage = 0.001,
             cv.folds = 5, n.cores = 6)
  best <- which.min(fit$cv.error)
  fit <- gbm(formula = Train_label[[classVar]] ~ .,
             data = as.data.frame(Train_set),
             distribution = 'bernoulli',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001, n.cores = 8)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunXGBoost <- function(Train_set, Train_label, mode, classVar){
  indexes = createFolds(Train_label[[classVar]], k = 5, list=T)
  CV <- unlist(lapply(indexes, function(pt){
    dtrain = xgb.DMatrix(data = Train_set[-pt, ],
                         label = Train_label[[classVar]][-pt])
    dtest = xgb.DMatrix(data = Train_set[pt, ],
                        label = Train_label[[classVar]][pt])
    watchlist <- list(train=dtrain, test=dtest)

    bst <- xgb.train(data=dtrain,
                     max.depth=2, eta=1, nthread = 2, nrounds=10,
                     watchlist=watchlist,
                     objective = "binary:logistic", verbose = F)
    which.min(bst$evaluation_log$test_logloss)
  }))

  nround <- as.numeric(names(which.max(table(CV))))
  fit <- xgboost(data = Train_set,
                 label = Train_label[[classVar]],
                 max.depth = 2, eta = 1, nthread = 2, nrounds = nround,
                 objective = "binary:logistic", verbose = F)
  fit$subFeature = colnames(Train_set)

  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunNaiveBayes <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  fit <- naiveBayes(eval(parse(text = paste(classVar, "~."))),
                    data = data)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# DRF不适用于二分类情况，因此删去
# RunDRF <- function(Train_set, Train_label, mode, classVar){
#   Train_label <- data.frame(
#     "0" = as.numeric(Train_label == 0),
#     "1" = as.numeric(Train_label == 1)
#   )
#   fit <- drf(X = Train_set,
#              Y = Train_label,
#              compute.variable.importance = F)
#   fit$subFeature = colnames(Train_set)
#
#   summary(predict(fit, functional = "mean", as.matrix(Train_set))$mean)
#
#   if (mode == "Model") return(fit)
#   if (mode == "Variable") return(ExtractVar(fit))
# }

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

standarize.fun <- function(indata, centerFlag, scaleFlag) {
  scale(indata, center=centerFlag, scale=scaleFlag)
}


#' Extract Variables Used in a Model
#'
#' @param fit A fitted model object. The function supports various model types such as `lognet`,
#' `glm`, `svm.formula`, `train`, `glmboost`, `plsRglmmodel`, `rfsrc`, `gbm`, `xgb.Booster`,
#' and `naiveBayes`.
#'
#' @return A character vector of selected feature names used by the model. For models that do not
#' perform variable selection, the function returns all variables. `(Intercept)` or `Intercept`
#' terms are excluded from the output.
#' @export
#'
#' @examples
#' # Example with a random forest model (rfsrc)
#' library(randomForestSRC)
#' fit <- rfsrc(Species ~ ., data = iris)
#' variables <- ExtractVar(fit)
#' print(variables)
#'
#' # Example with a generalized linear model (glm)
#' fit_glm <- glm(Species ~ Sepal.Length + Sepal.Width, data = iris, family = binomial)
#' variables_glm <- ExtractVar(fit_glm)
#' print(variables_glm)
ExtractVar <- function(fit){
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet" = rownames(coef(fit))[which(coef(fit)[, 1]!=0)], # 本身没有筛选变量的功能，但是可以舍去模型中系数为0的变量
    "glm" = names(coef(fit)), # 逐步回归可以对变量进行筛选
    "svm.formula" = fit$subFeature, # SVM对变量没有筛选功能，所以默认使用所有变量
    "train" = fit$coefnames, # LDA不能筛选变量，所以默认使用所有变量
    "glmboost" = names(coef(fit)[abs(coef(fit))>0]), # glmboost同样不具备筛选变量的功能，因此舍去模型中系数为0的变量
    "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs!=0], # plsRglmmodel同样不具备筛选变量的功能，因此舍去模型中系数为0的变量
    "rfsrc" = var.select(fit, verbose = F)$topvars, # rfsrc可以对变量进行筛选
    "gbm" = rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0], # gbm通过舍去重要性为0的变量来进行变量筛选
    "xgb.Booster" = fit$subFeature, # xgboost没有筛选变量的能力， 默认使用所有变量
    "naiveBayes" = fit$subFeature # naiveBayes没有筛选变量的能力，默认使用所有变量
    # "drf" = fit$subFeature # drf自带的变量系数提取函数会输出NA，因此默认使用所有变量
  ))

  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  return(Feature)
}











