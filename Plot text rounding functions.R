# Functions for rounding model text to 2 d.p. for putting on plots
# Written by Esther Dale, Feb 2021

#calculates p-value from model object f stat, rounds
p.value <- function(x){
  y <- summary(x)
  p.v <- pf(y$fstatistic[1], y$fstatistic[2], y$fstatistic[3], lower.tail = F)
  p.v <- round(p.v, digits = 2)
  if(nchar(p.v)==3){
    p.v <- paste(p.v, "0", sep="")
  }
  return(p.v)
}

# rounds R squared
r.sqr <- function(x){
  y <- summary(x)
  r.value <- round(y$r.squared, digits = 2)
  if(nchar(r.value)==3){
    r.value <- paste(r.value, "0", sep="")
  }
  return(r.value)
}

#rounds and combines df and f for putting on a plot
f.fig <- function(x){
  y <- summary(x)
  f <- round(y$fstatistic[1], digits=2)
  if(nchar(f)==3){
    f <- paste(f, "0", sep="")
  }
  
  df1 <- round(y$fstatistic[2], digits = 0)
  df2 <- round(y$fstatistic[3], digits = 0)
  f.text <- paste("F(", df1,",", df2, ") = ", f, sep="")
  return(f.text)
}
