# scatterWithHists <- function(x, y, histCols=c("lightblue","lightblue"), lhist=20, xlim=range(x), ylim=range(y), ...){
#   ## set up layout and graphical parameters
#   layMat <- matrix(c(1,4,3,2), ncol=2)
#   layout(layMat, widths=c(5/7, 2/7), heights=c(2/7, 5/7))
#   ospc <- 0.5                                                  # outer space
#   pext <- 4                                                    # par extension down and to the left
#   bspc <- 1                                                    # space between scatter plot and bar plots
#   par. <- par(mar=c(pext, pext, bspc, bspc), oma=rep(ospc, 4)) # plot parameters
#
#   ## barplot and line for x (top)
#   xhist <- hist(x, breaks=seq(xlim[1], xlim[2], length.out=lhist), plot=FALSE)
#   par(mar=c(0, pext, 0, 0))
#   barplot(xhist$density, axes=FALSE, ylim=c(0, max(xhist$density)), space=0, col=histCols[1])
#
#   ## barplot and line for y (right)
#   yhist <- hist(y, breaks=seq(ylim[1], ylim[2], length.out=lhist), plot=FALSE)
#   par(mar=c(pext, 0, 0, 0))
#   barplot(yhist$density, axes=FALSE, xlim=c(0, max(yhist$density)), space=0, col=histCols[2], horiz=TRUE)
#
#   ## overlap
#   dx <- density(x)
#   dy <- density(y)
#   par(mar=c(0, 0, 0, 0))
#   plot(dx, col=histCols[1], xlim=range(c(dx$x, dy$x)), ylim=range(c(dx$y, dy$y)),
#        lwd=4, type="l", main="", xlab="", ylab="", yaxt="n", xaxt="n", bty="n"
#   )
#   points(dy, col=histCols[2], type="l", lwd=3)
#
#   ## scatter plot
#   par(mar=c(pext, pext, 0, 0))
#   plot(x, y, xlim=xlim, ylim=ylim, ...)
# }
# eval(parse(text=factor))


