#calculates AIC

Aic <- function(x,y){
  aic <- -(2*x)+(2*y)
  return(aic)
}