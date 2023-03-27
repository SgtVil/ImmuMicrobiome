#' @title Collapse IgASeq
#'
#' @description Collapse the object returned by slide_z(). The slide_z() function needs to be performed in a for loop or lapply and be collapsed into a single object.
#'
#'
#' @param IgAseq Object returned by the slide_z() function
#'
#' @return IgAseq object. Refers to [slide_z] for more information on the IgAseq object.
#' @export
#'
#' @examples URL to come.
collapse_IgAseq = function(IgAseq){

  all= list()
  up= list()
  down= list()
  for(i in names(IgAseq)){
    if(class(IgAseq[[i]])=="IgASeq"){
      x1=IgAseq[[i]]
      all[[i]]= slot(x1, "ig_seq_all")

      up[[i]] = slot(x1, "ig_up")

      down[[i]] = slot(x1, "ig_down")
    }
  }

  # all= lapply(names(IgAseq), FUN= function(x)slot(IgAseq[[x]], "ig_seq_all"))
  all= do.call("rbind", all)
  # up = lapply(names(IgAseq), FUN= function(x)slot(IgAseq[[x]], "ig_up"))
  up= do.call("rbind", up)
  # down= lapply(names(IgAseq), FUN= function(x)slot(IgAseq[[x]], "ig_down"))
  down= do.call("rbind", down)
  res <- ig_seq(ig_seq_all = all,
                ig_up = up,
                ig_down = down)
  res
}
