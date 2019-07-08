#' a function that calculate the column cross production terms for a matrix or dataframe
#'
#' @param mydata 
#'
#' @return the column cross production terms from \code{\link{return}}
#' @export
#' @examples
#' x <- data.frame(runif(20,1,5),runif(20,1,5),runif(20,1,5))
#' rowprod(x)
rowprod <- function(mydata){
  
  re <- NULL
  
  for(i in 1 : ncol(mydata)){
    
    for(j in 1 : ncol(mydata)){
      
      re <- cbind(re, 0.5 * mydata[, i] * mydata[, j])
      
    }
  }
  return(re)
}

#' a function that calculate the column cross production terms for 2 matrices or dataframes
#'
#' @param mydata1 matrix or dataframe 1
#' @param mydata2 matrix or dataframe 2
#'
#' @return the column cross production terms from \code{\link{return}}
#' @export
#'
#' @examples
#' x <- data.frame(runif(20,1,5),runif(20,1,5),runif(20,1,5))
#' y <- data.frame(runif(20,1,5),runif(20,1,5),runif(20,1,5))
#' rowprod2(x,y)
rowprod2 <- function(mydata1, mydata2){
  
  re <- NULL
  
  for(i in 1 : ncol(mydata1)){
    
    for(j in 1 : ncol(mydata2)){
      
      re <- cbind(re, mydata1[,i] * mydata2[,j])
      
    }
  }
  return(re)
}

#' calculate the 1st order derivative of cross products
#'
#' @param x 
#'
#' @return the 1st order production terms from \code{\link{return}}
#' @export
#'
#' @examples \dontrun{
#' # No need to run this
#' }
#' x <- data.frame(runif(20,1,5),runif(20,1,5),runif(20,1,5))
#' crossderiv1(x)
crossderiv1 <- function(x){
  
  if(is.matrix(x)){x <- as.data.frame(x)}
  
  diagx <- diag(1,ncol(x))
  
  consx <- NULL
  
  for (m in 1 : ncol(x)){
    
    consxk <- NULL
    
    for (k in 1: nrow(x)){
      
      consxk.temp <- as.numeric(0.5 * (t(x[k,]) %*% diagx[m,] + t(t(x[k,]) %*% diagx[m,])))
      
      consxk <- rbind(consxk, consxk.temp)
    }
    
    consx <- rbind(consx, consxk)
    
    
  }
  
  return(consx)
  
}

#' calculate the 1st order derivative of cross products of 2 matrices or dataframes
#'
#' @param x the first matrix
#' @param y the second matrix
#'
#' @return the 1st order production terms from \code{\link{return}}
#' @export
#'
#' @examples \dontrun{
#' # No need to run this
#' }
#' x <- data.frame(runif(20,1,5),runif(20,1,5),runif(20,1,5))
#' y <- data.frame(runif(20,1,5),runif(20,1,5))
#' crossderiv1(x)
crossderiv2 <- function(x,y){
  
  if(!is.data.frame(x) || !is.data.frame(y)){
    
    x <- as.data.frame(x)
    
    y <- as.data.frame(y)
    
  }
  if(nrow(x) != nrow(y))
    stop("the 2 data sets should have the same number of rows")
  
  consxy <- NULL
  
  for (j in 1: ncol(x)){
    
    cons <- NULL
    
    diagx <- diag(1,ncol(x))
    
    for (i in 1:ncol(y)){
      
      cons.temp <- y[,i] %*% t(diagx[,j])
      
      cons <- cbind(cons,cons.temp)
      
    }
    
    consxy <- rbind(consxy, cons)
    
  }
  
  return(consxy)
  
}

#' The function to derive the symmetric constraints
#'
#' @param n 
#'
#' @return the symmetry constraints from \code{\link{return}}
#' @export
#'
#' @examples \dontrun{
#' # No need to run this
#' }
#' conssym(3)
conssym <- function(n){
  
  if(n > 1){
    
    conssym <- NULL
    
    for (i in 1:(n-1)) {
      
      for (j in (i+1):n) {
        
        matrix.temp <- matrix(0,n,n)
        
        matrix.temp[i,j] <- 1
        
        matrix.temp[j,i] <- -1
        
        numeric.matrix.temp <- as.numeric(t(matrix.temp))
        
        conssym <- rbind(conssym, numeric.matrix.temp)
        
      }
      
    }
    
    return(conssym)
    
  }else{return(matrix(NA,0,0))}
}

#' A function that separates the target function and the constraints
#'
#' @param mydata A matrix
#'
#' @return the separated matrix from \code{\link{return}}
#' @export
#'
#' @examples  \dontrun{
#' # No need to run this
#' }
#' x <- cbind(runif(20,1,5),runif(20,1,5),runif(20,1,5))
#' change2(x)
change2 <- function(mydata){
  
  J <- ncol(mydata)
  
  z <- matrix(NA, nrow(mydata),ncol = 2*J)
  
  for(j in 1 : J){
    
    z[, 2*j - 1] <- mydata[,j]
    
    z[, 2*j] <- -mydata[,j]
    
  }
  return(z)
}


