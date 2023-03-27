#' Create a S4 object for the IgASeq data
#'
#' @slot ig_seq_all data.frame.
#' @slot ig_up data.frame.
#' @slot ig_down data.frame.
#' @slot ellipse_data data.frame.
#' @slot taxa data.frame.
#'
#' @return A class IgASeq
#' @export
#'
#' @examples .
ig_seq <- setClass(Class = "IgASeq",package = "ImmuMicrobiome",
                   slots = list(ig_seq_all = "data.frame",
                                ig_up = "data.frame",
                                ig_down = "data.frame",
                                ellipse_data= "data.frame",
                                taxa= "data.frame"))
