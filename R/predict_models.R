#' Predict the classification based on the model training.
#' @param model_obj A object returned by [screen_model()].
#'
#' @return A list of predicted values and the confusion matrix for each variables.
#' @export
#'
#' @examples TBD
predict_models = function(model_obj){
  if(is.null(model_obj$pred_y)){
    stop("your data wasn't trained, use train=T in predict_models")
  } else {
    res= mclapply(1:length(model_obj), function(x){
      predicted = model_obj$model[[x]]
      true_data = model_obj$pred_x[[x]]
      res= stats::predict(predicted, true_data)
    })

    confusion = lapply(1:length(model_obj), function(x){
      predicted = res[[x]]
      true_data = model_obj$pred_y[[x]]
      confusionMatrix(predicted, factor(true_data))
    })
  }

  return(list(predicted= res, conf_matrix = confusion))
}
