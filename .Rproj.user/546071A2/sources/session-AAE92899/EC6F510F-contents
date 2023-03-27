#' Create a S4 object for the machine learning results.
#'
#' @slot model List of models returned
#' @slot pred List of predictions
#' @slot conf List of confusion matrix and the attached statistics
#' @slot data List of the data used to generate this object
#'
#' @return A class ml
#' @export
#'
#' @examples .
ml <- setClass(Class = "ml", package = "ImmuMicrobiome",
                   slots = c("res_model", "bacterial_loading", "res_prediction", "res_conf_matrix", "data"))
