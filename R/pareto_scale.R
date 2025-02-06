# pareto_scale = function(x){
#   rowmean = apply(x, 1, mean)
#   rowsd = apply(x, 1, sd)
#   rowsqrt = sqrt(rowsd)
#   rv = sweep(x, 1, rowmean, "-")
#   rv = sweep(rv, 1, rowsqrt, "/")
#   return(rv)
# }

# doi:10.1186/1471-2164-7-142
pareto_scale =  function(x){
  (x - mean(x)) / sqrt(sd(x))
}

cubic_root = function(x) {
  ifelse(x >= 0, x^(1/3), -((-x)^(1/3)))
}

auto_scale = function(x){
  (x - mean(x)) / sd(x)
}

level_scaling = function(x){
  (x - mean(x)) / mean(x)
}

vast_scaling= function(x){
  (x - mean(x)) / sd(x) * mean(x) / sd(x)
}

median_centering = function(x){
  x-median(x)
}

geometric_mean = function(x) x- exp(mean(log(x[x!=0])))

power_transf = function(x){
  s
}

sum_scaling = function(x)x/sum(x)
