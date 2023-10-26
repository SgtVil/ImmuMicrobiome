Shannon= function(physeq, abs=TRUE){
  if(abs==T){
    shannon <- function(x){
      x <- x[!is.na(x)]
      result <- 0
      nval <- length(x)
      sumval <- sum(x)

        if (sumval != 0) {
          for (i in 1:nval) {
            if (x[i] != 0) {
              result <- result + x[i]*log(x[i])
            }
          }
          result <- (log(sumval)-(1/sumval)*result)
        } else {
          return <- -1 #Error value
        }

    }
   as.data.frame( cbind(sample_data(physeq),apply(otu_table(physeq), 2, function(x)shannon(x)) ))
  }

   else{
     Shannon <- function(x){
       x <- x[!is.na(x)]
       result <- 0
       nval <- length(x)
       sumval <- sum(x)

       if (sumval != 0) {
         for (i in 1:nval) {
           if (x[i] != 0) {
             result <- result + x[i]*log(x[i])
           }
         }
         result <- (log(sumval)-(1/sumval)*result)/log(nval)
       } else {
         return <- -1 #Error value
       }

     }

    as.data.frame(cbind(sample_data, apply(otu_table(physeq), 2, function(x)shannon(x))))
   }
}
