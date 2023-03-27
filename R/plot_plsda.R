plot_model = function(obj, comp_x="comp1", comp_y="comp2", plot_name= "model.pdf", cores=2){

  var = sapply(obj$model, "[", "variates")
  X= sapply(var, "[", "X")
  Y = sapply(obj$model, "[", "Y")
  contrib = sapply(obj$model, "[", "prop_expl_var")

  names(X)= names(obj$model) ;  names(Y)= names(obj$model) ; names(contrib)= names(obj$model)
print(dim(contrib))

  plot_p = function(name){
    if(length(unique(Y[[name]]))>10){
      warning(paste(name, "wasn't plotted because there is too much groups"))
    } else {
      ggplot(data = as.data.frame(X[[name]]), aes_string(comp_x, comp_y, fill=Y[[name]]))+
        geom_point(shape=21)+
        stat_ellipse(data= as.data.frame(X[[name]]),
                     mapping = aes_string(comp_x, comp_y, fill=Y[[name]], color= Y[[name]]), geom="polygon", alpha=0.3, level = 0.9)+
        labs(x= paste("Component 1:", round(contrib[[name]][[1]][1]*100,2), "%"),
             y=paste("Component 2:", round(contrib[[name]][[1]][2]*100,2), "%"),
             color=name, fill=name)+
        theme(panel.grid = element_blank(),
              legend.position = "bottom")
    }

  }
 pdf(plot_name, width = 12, height = 10)
  mclapply(names(X), plot_p, mc.cores = cores)
dev.off()

  }



