#' Cross validate the training data.
#'
#' @param obj An object of returned by [screen_models()]
#' @param y Character vector of the desire variable to test.
#' @param model Type of model to use
#' @param repetition The number of repetition
#' @param cores The number of cores to use
#' @param ncomp The number of components you want to keep in plsda or splsda
#' @param repeats For repeated k-fold cross-validation only: the number of complete sets of folds to compute. Default = 3
#' @param prop The proportion of the variables to use in the training model. Default = 0.5.
#' @param number Either the number of folds or number of resampling iterations
#' @param ...
#'
#'
#' @import foreach
#' @return A list containing of ml class object.
#' @item{res_model}{The [caret::train()] value}
#' @item{res_prediction}{The prediction of the model}
#' @item{res_conf_matrix}{The confusion matrix and associated values, see [caret::confusionMatrix()].
@item{data}{The data used in each models}
#'
#' @export
#' @import doParallel
#' @import foreach
#' @examples TBD
cross_validate_data = function(obj, y, model="plsda", repetition= 1000, cores=1, ncomp=10, number=10,
                               repeats=3,  prop= 0.5, ...){

  data= obj[[y]]@data$X
  variables= obj[[y]]@data$Y

  control <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 3, search = "random")
 mclapply(1:repetition,
                       function(x){

                         index = caret::createDataPartition(variables, list = F, p = prop)
                         train=NULL
                         train$y= variables[index]
                         train$x= data[index,]

                         pred= NULL
                         pred$y= variables[-index]
                         pred$x= data[-index,]

                         res_model=NULL

                         if(model=="plsda")  {
                           res_model = caret::plsda(x = train$x,  y= as.factor( train$y), ncomp=20, ...)
                           loads = varImp(res_model)
                         } else  if(model=="glmnet"){
                           res_model= caret::train(train$x, train$y, method="glmnet", trControl=control, ...)
                           loads = varImp(res_model)[["importance"]]
                         } else if(model=="rf") {
                           message("Random Forest algorithm takes to long to cross validate like this")
                         }

                         pred_model= NULL
                         conf= NULL

                         pred_model= predict(res_model, pred$x)
                         conf= caret::confusionMatrix(pred_model, factor(pred$y))

                       res = ml(res_model =  res_model,
                                bacterial_loading = loads,
                                res_prediction =  pred_model,
                                res_conf_matrix = conf,
                                data = list(data=NULL,
                                            train= train,
                                            pred= pred))

                       }, mc.cores = cores )
  }
