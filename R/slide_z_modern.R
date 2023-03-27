#' Modern version of the Z test on the IgASeq data.
#' @description \code{slide_z_modern} will use the dispersion between the two negative fractions to mesure the technical dispersion and apply the Z stat on the pos/neg1 ration
#'
#' @param df The dataframe returned by \code{\link{seq_table}}
#' @param log2_ratio The \code{\link{log_ratio}} of neg
#' @param log10_abundance The
#' @param deltaX The
#' @param alpha The
#'
#' @return A slidez
#' @export
#'
#' @examples No example
slide_z_modern <- function(df, log2_ratio, log10_abundance, deltaX=deltaX, alpha=alpha){


  backup <- as.data.frame(df)

  # if(dim(backup)[2]<5 ) warning(paste(colnames(df),"There is multiple occurrence of positive and negative references in your data"))
  'Calculate sliding Z-score with overlap'
  backup2 <- backup
  Overlap <- deltaX/2
  Cycles <- nrow(backup2) %/% deltaX  #(Integer division)
  Rest <- nrow(backup2) %% deltaX

  if (length(backup[,1])!=0) {

    #(modulus)
    if ((Rest>Overlap) | (Cycles==0)){
      Cycles <- Cycles+1
    }
    for (j in 1:Cycles) {
      if (j == 1 && j!=Cycles) {
        a <- 1
        b <- deltaX+Overlap
      } else if (j>1 && j<Cycles) {
        a <- (j-1)*deltaX+1-Overlap
        b <- j*deltaX+Overlap
      } else if(j>1 && j==Cycles){
        a <- (j-1)*deltaX+1-Overlap
        b <- nrow(backup2)
      } else if(j==1 && j==Cycles){
        a <- 1
        b <- nrow(backup2)
      }
      DataTemp <- backup2[a:b,]

      # Centre reduction

      SlideNormTemp <- DataTemp$log2_ratio - mean(DataTemp$log2_neg_ratio)
      SlideZTemp <- (DataTemp$log2_ratio - mean(DataTemp$log2_neg_ratio))/ sd(DataTemp$log2_neg_ratio - mean(DataTemp$log2_neg_ratio))
      # Merge analyzed data
      if (j == 1 && j!=Cycles) {
        a <- 1
        b <- deltaX
        SlideNormTemp <- SlideNormTemp[a:b]
        SlideZTemp <- SlideZTemp[a:b]
        SlideNormFinal <- SlideNormTemp
        SlideZFinal <- SlideZTemp
      } else if (j>1 && j<Cycles){
        a <- Overlap+1
        b <- Overlap+deltaX
        SlideNormTemp <- SlideNormTemp[a:b]
        SlideZTemp <- SlideZTemp[a:b]
        SlideNormFinal <- c(SlideNormFinal, SlideNormTemp)
        SlideZFinal <- c(SlideZFinal, SlideZTemp)
      } else if(j>1 && j==Cycles) {
        a <- Overlap+1
        b <- length(SlideNormTemp)
        SlideNormTemp <- SlideNormTemp[a:b]
        SlideZTemp <- SlideZTemp[a:b]
        SlideNormFinal <- c(SlideNormFinal, SlideNormTemp)
        SlideZFinal <- c(SlideZFinal, SlideZTemp)
      } else if(j==1 && j==Cycles){
        a <- 1
        b <- length(SlideNormTemp)
        SlideNormTemp <- SlideNormTemp[a:b]
        SlideZTemp <- SlideZTemp[a:b]
        SlideNormFinal <- SlideNormTemp
        SlideZFinal <- SlideZTemp
      }
    }
    backup2$SlideNorm <- SlideNormFinal
    backup2$score <- SlideZFinal

    Zalpha <- qnorm(p=(1-alpha/2),mean=0, sd=1)
    ig_up <- backup2[with(backup2,which(score > Zalpha)), ]
    ig_down <- backup2[with(backup2,which(score < -Zalpha)), ]
    res <- list(ig_seq_all=backup2, ig_up=ig_up, ig_down=ig_down)

  }
}
