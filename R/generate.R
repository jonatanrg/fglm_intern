generate_normal <- function(X, b, d = 1, ymax = 1){
  # Do argument checking
  stopifnot(is.matrix(X))
  p <- ncol(X)
  n <- nrow(X)
  stopifnot(is.numeric(b), length(b) == p)
  stopifnot(is.numeric(d), length(d) == 1, d > 0)
  stopifnot(is.numeric(ymax), length(ymax) == 1, ymax >= 0)
  
  y <- rnorm(n = nrow(X),0,1)
  for (i in 1:length(y)){
    if (y[i] > -X[i,]%*%b){
      y[i] <- 1
    }
    else{
      y[i] <- 0
    }
  }
  return(y)
}