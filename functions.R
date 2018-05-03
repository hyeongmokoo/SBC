SBC <- function(x.e, x.v, nb)
{
  # purpose: Calculating spatial autocorrelation incorporating uncertainty information
  # Arguments:
  #   x.e: a vector of estimates
  #   x.v: a vector of uncertainty values, forms of standard errors
  #   nb: nb object
  # Fuction called:
  #   Bhatta.dist, Bhatta.coef, He.dist
  # Author: Hyeongmo Koo (05/03/2018)
  
  
  n <- length(x.e)
  
  res.bd <- vector(mode = 'numeric', length = n)
  res.bc <- vector(mode = 'numeric', length = n)
  res.hd <- vector(mode = 'numeric', length = n)
  
  nb.n <- card(nb)
  
  for(i in 1:n){
    tmp.bd <- 0
    tmp.bc <- 0
    tmp.hd <- 0
    
    for (j in 1:nb.n[i])
    {
      nb.idx <- nb[[i]][j]
      lc.bd <- Bhatta.dist(x.e[i], x.e[nb.idx], x.v[i], x.v[nb.idx])
      lc.bc <- Bhatta.coef(x.e[i], x.e[nb.idx], x.v[i], x.v[nb.idx])
      lc.hd <- He.dist(x.e[i], x.e[nb.idx], x.v[i], x.v[nb.idx])
      
      tmp.bd <- tmp.bd + lc.bd
      tmp.bc <- tmp.bc + lc.bc
      tmp.hd <- tmp.hd + lc.hd
    }
    
    res.bd[i]<-tmp.bd/nb.n[i]
    res.bc[i]<-tmp.bc/nb.n[i]
    res.hd[i]<-tmp.hd/nb.n[i]
  }
  global.res <- c(mean(res.bd), mean(res.bc), mean(res.hd))
  names(global.res) <- c("bd", "bc", "hd")
  res <- list(global.res, res.bd, res.bc, res.hd)
  names(res) <- c("global", "local.bd", "local.bc", "local.hd")
  res
}



Bhatta.dist <- function(x1, x2, v1, v2){
  
  var.comp <- log((1/4)*((v1^2/v2^2)+(v2^2/v1^2)+2))
  mean.var.comp <- ((x1-x2)^2)/(v1^2+v2^2)
  bha.dist <- (1/4)*(var.comp+mean.var.comp)
  
  return(bha.dist)
}
Bhatta.coef <- function(x1, x2, v1, v2){
  
  var.comp <- log((1/4)*((v1^2/v2^2)+(v2^2/v1^2)+2))
  mean.var.comp <- ((x1-x2)^2)/(v1^2+v2^2)
  bha.dist <- (1/4)*(var.comp+mean.var.comp)
  bha.coef <- exp(bha.dist * -1)
  return(bha.coef)
}

He.dist <- function(x1, x2, v1, v2){
  
  var.comp <- log((1/4)*((v1^2/v2^2)+(v2^2/v1^2)+2))
  mean.var.comp <- ((x1-x2)^2)/(v1^2+v2^2)
  bha.dist <- (1/4)*(var.comp+mean.var.comp)
  bha.coef <- exp(bha.dist * -1)
  hell.dist <- sqrt(1-bha.coef)
  return(hell.dist)
}

