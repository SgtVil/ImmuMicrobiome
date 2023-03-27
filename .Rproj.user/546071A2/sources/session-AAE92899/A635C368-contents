#' @export
prepare_physeq = function(obj,  prop=0.5, train=T, ... ){
  # extract phyloseq sub objects
  if(taxa_are_rows(obj)){
    obj = reverseASV(obj)
  }
  x = as(otu_table(obj), "matrix")
  x= x[, apply(otu_table(obj), 2, function(x)sum(x)>0)]
  y= as(phyloseq::sample_data(obj), "data.frame")

  data = pre_process(x, y)

  train_y=NULL
  train_x=NULL
  pred_y = NULL
  pred_x = NULL
  # prepare training and testing data
  if(train== T){
    for(i in names(data$y)){
      index = caret::createDataPartition(data$y[[i]], list = F, p = prop)
      train_y[[i]]= data$y[[i]][index]
      train_x[[i]]= data$x[[i]][index,]

      pred_y[[i]]= data$y[[i]][-index]
      pred_x[[i]]= data$x[[i]][-index,]
    }
  }
  if(train==F){
    for(i in names(data$y)){
      train_y[[i]]= data$y[[i]]
      train_x[[i]]= data$x[[i]]
    }
  }
return(list(train_y= train_y, train_x= train_x, pred_x =pred_x, pred_y= pred_y))
}
