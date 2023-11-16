#' Title Screen a model on all your variables.
#' @description A function to test machine learning on your different variables of the metadata.
#' This function supports phyloseq object as input or a list containing a values matrix and a variable dataframe.
#' This function is based on [caret].
#' @param obj An object of class phyloseq or a list of two levels containing a matrix of values and a dataframe of variables.
#' @param model The model tested, for now this function supports [randomForest::randomForest()], [glmnet::glmnet()] and [caret::plsda()].
#' @param cores The number of cores you want to use.
#' @param ncomp Number of component you want to keep for [caret::plsda()]. Default =10.
#'
#' @param number Either the number of folds or number of resampling iterations
#' @param repeats For repeated k-fold cross-validation only: the number of complete sets of folds to compute. Default = 3
#' @param prop The proportion of the variables to use in the training model. Default = 0.5.
#' @param train Wheter or not you want to train the models in order to predict.
#' If you select FALSE the function will use all the data to perform classification or find
#' the best features that explain your data according to the variables. Default=T.
#' @param ... Argument passed to [caret::train()].
#'
#' @return A list containing
#' \item{model}{A [caret::train()] class}
#' \item{Y}{The trimmed variables}
#' \item{X}{The trimmed values}
#' \item{train_y}{The variables used to train the model}
#' \item{train_x}{The values used to train the model}
#' \item{pred_x}{The predicted values}
#' \item{pred_y}{The predicted classification}
#'
#'
#' @examples
#' screen_models(enterotype, model="glmnet", prop=0.5, ncomp=20)
#'
#' screen_models(enterotype, model="plsda", cores=1, ncomp=10, number=10, repats=5, prop=0.6, train=T)
screen_models = function(obj, model="plsda", cores=1, ncomp=10, number=10,
                         repeats=3,  prop= 0.5, train= T, ... ){

  if(class(obj)=="phyloseq"){
    prepare_physeq(obj, model="plsda", ncomp=10, prop= 0.5, train= T, ...)
  } else if(class(obj)=="matrix" | class(obj)=="data.frame"){
    prepare_non_physeq(obj, model="plsda", ncomp=10, prop= 0.5, train= T, ...)
  }
  # pred_model= NULL
  # conf= NULL

  if(model == "plsda"){ # run plsda
    res_model= mclapply(X = 1:length(train_y),
                  FUN = function(i){
                    x = train_x[[i]]
                    y = train_y[[i]]

                    res = caret::plsda(x = train_x[[i]],  y= as.factor(train_y[[i]]), ncomp=ncomp)

                    pred_model= predict(res, pred_x[[i]])
                    conf= caret::confusionMatrix(pred_model, factor(pred_y[[i]]))
                    loads = varImp(res)

                    res_model = ml(res_model = res,
                                   bacterial_loading = loads,
                                   res_prediction = pred_model,
                                   res_conf_matrix = conf,
                                   data= list(X= x,
                                              Y= y,
                                              train_y= train_y[[i]],
                                              train_x= train_x[[i]],
                                              pred_x= pred_x[[i]],
                                              pred_y= pred_y[[i]]))
                    },
                  mc.cores = cores)

  } else if(model=="splsda"){
    res_model= mclapply(X = 1:length(train_y),
                  FUN = function(i){
                    x = train_x[[i]]
                    y = train_y[[i]]

                   res= caret::splsda(x = train_x[[i]],
                                  y= as.factor(train_y[[i]]),
                                  ncomp=ncomp, ...)
                    pred_model= predict(res, pred_x[[i]])
                    conf= caret::confusionMatrix(pred_model, factor(pred_y[[i]]))
                    loads = varImp(res)

                    res_model = ml(res_model = res,
                                   bacterial_loading = loads,
                                   res_prediction = pred_model,
                                   res_conf_matrix = conf,
                                   data= list(X= x,
                                              Y= y,
                                              train_y= train_y[[i]],
                                              train_x= train_x[[i]],
                                              pred_x= pred_x[[i]],
                                              pred_y= pred_y[[i]]))
                    },
                  mc.cores = cores)

  }  else {
    control = trainControl(method="repeatedcv", number = number, repeats=repeats)
    res_model = mclapply(X= 1:length(train_y),
                   FUN = function(i){

                    res= caret::train(x=train_x[[i]], y=train_y[[i]], method= model, metric="Accuracy")
                    pred_model= predict(res, pred_x[[i]])
                    conf= caret::confusionMatrix(pred_model, factor(pred_y[[i]]))
                    loads = varImp(res)[["importance"]]

                    res_model = ml(res_model = res,
                                   bacterial_loading = loads,
                                   res_prediction = pred_model,
                                   res_conf_matrix = conf,
                                   data= list(X= x,
                                              Y= y,
                                              train_y= train_y[[i]],
                                              train_x= train_x[[i]],
                                              pred_x= pred_x[[i]],
                                              pred_y= pred_y[[i]]))
                    },
                   mc.cores= cores)
  }

  names(res_model)= names(data$y)
 return(res_model)
}
