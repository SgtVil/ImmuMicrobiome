#' @export

pre_process = function(x, y){

    bool = sapply(y, unique)

  new_y = NULL
  new_x = NULL

  for(col in 1:length(colnames(y))){ #check the variables before running plsda
    # print(col)
    if(length(bool[[col]])==1 | ((anyNA(y[,col]) | any(str_detect(y[,col], "NA"))) & length(bool[[col]])==2)){ # is the col as only one condition
      print(paste(colnames(y)[col], ": vector of 1 condition"))

    } else if(anyNA(y[,col]) & length(bool[[col]])>=2){ # remove NA
      # print(paste(colnames(y)[col], "contains NA"))
      new_x[[colnames(y)[col]]] = x[!is.na(y[,col]),]
      new_y[[colnames(y)[col]]] = y[!is.na(y[,col]), col]

    } else if(any(str_detect(y[,col], "NA")) & length(bool[[col]])>=2){ # remove "NA" values
      new_x[[colnames(y)[col]]] = x[!str_detect(y[,col], "NA"),]
      new_y[[colnames(y)[col]]] = y[!str_detect(y[,col], "NA"), col]

    } else if(length(bool[[col]])==length(y[,col])){
      print(paste(colnames(y)[col], ": is constituted of unique observations"))

    } else if (length(bool[[col]])>=2 & !anyNA(y[,col])) {
      # print(paste(colnames(y)[col],"is good to go"))
      new_y[[colnames(y)[col]]] = y[,col] # make use only of the columns that a good for plsda
      new_x[[colnames(y)[col]]] = x
    }


  }
  return(list(x= new_x, y= new_y))
}
