pareto_scale = function(x){
  rowmean = apply(x, 1, mean)
  rowsd = apply(x, 1, sd)
  rowsqrt = sqrt(rowsd)
  rv = sweep(x, 1, rowmean, "-")
  rv = sweep(rv, 1, rowsqrt, "/")
  return(rv)
}
