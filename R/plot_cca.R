plot_cca <- function(cca.res, axis_x=1, axis_y=2, color_vector= c("cyan4","brown","deepskyblue", "black","red"),
                     col_vector, draw= "lines",
                     ylimits="auto", xlimits= "auto", lwd=1, conf=0.9, cex=2, text=F, ncol=1,
                     x.intersp = 1, y.intersp=0.5, where="topleft", inset=0.2, stat.cex=2){
  old.par = par()
  on.exit(layout(matrix(c(1,1))))
  
  
  p_li= scores(cca.res, display = "sites", choices=c(axis_x, axis_y) )
  
  ca1 = scores(cca.res, display = "species", scaling="species")[,1]
  ca2 = scores(cca.res, display = "species", scaling="species")[,2]
  
  
  if(xlimits=="auto"){
    xlimits = c(p_li[!is.na(p_li[,axis_x]),axis_x] %>% max  /0.5, p_li[!is.na(p_li[,axis_x]),axis_x] %>% min  /0.8)
  }
  
  if(ylimits=="auto"){
    ylimits = c(p_li[!is.na(p_li[,axis_y]),axis_y]  %>% max  /0.5, p_li[!is.na(p_li[,axis_y]),axis_y]  %>% min  /0.8)
  } else{
    ylimits= ylim
    xlimits= ylim
  }
  # prepare the matrix
  layout(matrix(c(2,2,2,4,
                  1,1,1,3,
                  1,1,1,3,
                  1,1,1,3),
                nrow = 4,
                ncol = 4,
                byrow = TRUE))
  
  par(mar=c(0,0,0,0))
  
  #1
  plot(p_li[,axis_x], p_li[,axis_y], bg= col2, axes=F, xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)
  disp= ordispider(p_li, groups = fac, col = adjustcolor(color_vector, alpha=0.3), lwd=lwd,  ylim=ylimits, xlim=xlimits)
  ordiellipse(p_li, groups= fac, conf= conf, col = color_vector, lwd = lwd, draw= draw,  ylim=ylimits, xlim=xlimits)
  points(p_li[,axis_x], p_li[,axis_y], bg= col2, xlab="", ylab="", las=2, pch=21, cex=cex, ylim=ylimits, xlim=xlimits, ...)
  
  if(text){
    text(x= unique(disp[,1:2]),  labels=unique(group), col="black", cex=cex, font=font)
  }
  
  if(stat=="permanova" | stat=="envfit"){
    
    legend(where, inset = 0.2, legend=p.val, bty="n", cex= stat.cex)
  }
  #2
  boxplot(p_li[,axis_x]~fac, data=p_li, horizontal=T, axes=F,  xlab=NULL, ylab=NULL,
          col= color_vector , xaxt="n", lwd=lwd/2,  ylim=xlimits)
  stripchart(p_li[,axis_x]~fac, data=p_li, method = "jitter", vertical=F, add=T, pch=21,
             col= "black", bg="gray", lwd=lwd/2,  ylim=xlimits)
  #3
  boxplot(p_li[,axis_y]~fac, data=p_li, axes=F,  ylab="", xlab="", col= color_vector, lwd=lwd/2,  ylim=ylimits)
  stripchart(p_li[,axis_y]~fac, data=p_li, method = "jitter", vertical=T, add=T, pch=21,
             col= "black", bg="gray", lwd=lwd/2,  ylim=ylimits)
  
  # #4
  plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="", axes=F)
  legend("center", legend=unique(fac), col = col1, title= legend_title, pch= 20,
         cex= cex, bty="n", ncol= ncol, y.intersp = y.intersp)
  
  # return(head(p$eig/sum(p$eig)*100,5))
  
}
