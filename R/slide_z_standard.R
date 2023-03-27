slide_z_standard <- function(df, log2_ratio, log10_abundance, positive_sorted_sample, negative_sorted_sample, deltaX=deltaX, alpha=alpha){

  backup <- as.data.frame(df)

  # if(dim(backup)[2]<5 ) warning(paste(colnames(df),"There is multiple occurrence of positive and negative references in your data"))
  'Calculate sliding Z-score with overlap'
  backup2 <- backup
  Overlap <- deltaX/2
  Cycles <- nrow(backup2) %/% deltaX  #(Integer division)
  Rest <- nrow(backup2) %% deltaX

  if (length(backup[,1])!=0) {
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
      SlideNormTemp <- DataTemp$log2_ratio - mean(DataTemp$log2_ratio)
      SlideZTemp <- (DataTemp$log2_ratio - mean(DataTemp$log2_ratio))/ sd(DataTemp$log2_ratio - mean(DataTemp$log2_ratio))
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

    # backup2$ASV <- row.names(backup2)
    Zalpha <- qnorm(p=(1-alpha/2),mean=0, sd=1)
    ig_up <- backup2[with(backup2,which(score > Zalpha)), ]
    ig_down <- backup2[with(backup2,which(score < -Zalpha)), ]
    # ig_up$ASV <- row.names(DiffUp)
    # ig_dow$ASV <- row.names(DiffDown)
    res <- list(ig_seq=backup2, ig_up=ig_up, ig_down=ig_down)

  }
}
