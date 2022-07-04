library(HRW)
library(mvtnorm)
library(Matrix)

# spline.des <- function (knots, x, ord = 4, derivs = integer(length(x)), outer.ok = FALSE, 
#           sparse = FALSE) 
# {
#   if (is.unsorted(knots <- as.numeric(knots))) 
#     knots <- sort.int(knots)
#   list(knots = knots, order = ord, derivs = derivs, design = splineDesign(knots, 
#                                                                           x, ord, derivs, outer.ok = outer.ok, sparse = sparse))
# }

# ZOSull1 <- function (x, range.x, intKnots, drv = 0) 
# {
#   if (drv > 2) 
#     stop("splines not smooth enough for more than 2 derivatives")
#   if (missing(range.x)) 
#     range.x <- c(1.05 * min(x) - 0.05 * max(x), 1.05 * max(x) - 
#                    0.05 * min(x))
#   if (missing(intKnots)) {
#     numIntKnots <- min(length(unique(x)), 35)
#     intKnots <- quantile(unique(x), seq(0, 1, length = (numIntKnots + 
#                                                           2))[-c(1, (numIntKnots + 2))])
#   }
#   numIntKnots <- length(intKnots)
#   allKnots <- c(rep(range.x[1], 4), intKnots, rep(range.x[2], 
#                                                   4))
#   K <- length(intKnots)
#   L <- 3 * (K + 8)
#   xtilde <- (rep(allKnots, each = 3)[-c(1, (L - 1), L)] + rep(allKnots, 
#                                                               each = 3)[-c(1, 2, L)])/2
#   wts <- rep(diff(allKnots), each = 3) * rep(c(1, 4, 1)/6, 
#                                              K + 7)
#   Bdd <- spline.des(allKnots, xtilde, derivs = rep(2, length(xtilde)), 
#                     outer.ok = TRUE)$design
#   Omega <- t(Bdd * wts) %*% Bdd
#   eigOmega <- eigen(Omega)
#   indsZ <- 1:(numIntKnots + 2)
#   UZ <- eigOmega$vectors[, indsZ]
#   LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
#   indsX <- (numIntKnots + 3):(numIntKnots + 4)
#   UX <- eigOmega$vectors[, indsX]
#   L <- cbind(UX, LZ)
#   # stabCheck <- t(crossprod(L, t(crossprod(L, Omega))))
#   # if (sum(stabCheck^2) > 1.0001 * (numIntKnots + 2)) 
#   #   print("WARNING: NUMERICAL INSTABILITY ARISING\\\n              FROM SPECTRAL DECOMPOSITION")
#   B <- spline.des(allKnots, x, derivs = rep(drv, length(x)), 
#                   outer.ok = TRUE)$design
#   Z <- B %*% LZ
#   attr(Z, "range.x") <- range.x
#   attr(Z, "intKnots") <- intKnots
#   return(Z)
# }

# O'sullivan
cal_z <- function(numBasis ,XsplPreds){
  
  numSplCompons <- ncol(XsplPreds)
  if (numSplCompons == 0) {
    Zspl_p <- NULL
    range.x.list <- NULL
    intKnots.list <- NULL
    ncZspl_p <- NULL
  }
  if (numSplCompons > 0) {
    range.x.list <- vector("list", numSplCompons)
    intKnots.list <- vector("list", numSplCompons)
    ncZspl_p <- NULL
    Zspl_p <- NULL
    for (jSpl in 1:numSplCompons) {
      xCurr <- XsplPreds[, jSpl]
      range.xVal <- c(min(xCurr), max(xCurr))
      numIntKnotsVal <- numBasis[jSpl] - 2
      intKnots <- as.numeric(quantile(unique(xCurr), seq(0, 1, length = numIntKnotsVal + 2)[-c(1, numIntKnotsVal + 2)]))
      Zcurr <- ZOSull(xCurr, intKnots = intKnots, range.x = range.xVal)
      range.x.list[[jSpl]] <- attr(Zcurr, "range.x")
      intKnots.list[[jSpl]] <- attr(Zcurr, "intKnots")
      ncZspl_p <- c(ncZspl_p, ncol(Zcurr))
      Zspl_p <- cbind(Zspl_p, Zcurr)
      
    }
  }
  return(Zspl_p)
}

cal_z_kp <- function(numBasis ,XsplPreds){
  
  numSplCompons <- ncol(XsplPreds)
  if (numSplCompons == 0) {
    Zspl_p <- NULL
    range.x.list <- NULL
    intKnots.list <- NULL
    ncZspl_p <- NULL
  }
  if (numSplCompons > 0) {
    range.x.list <- vector("list", numSplCompons)
    intKnots.list <- vector("list", numSplCompons)
    ncZspl_p <- NULL
    Zspl_p <- NULL
    for (jSpl in 1:numSplCompons) {
      xCurr <- XsplPreds[, jSpl]
      range.xVal <- c(min(xCurr), max(xCurr))
      numIntKnotsVal <- numBasis[jSpl] - 2
      intKnots <- as.numeric(quantile(unique(xCurr), 
                                      seq(0, 1, 
                                          length = numIntKnotsVal + 2)[-c(1, numIntKnotsVal + 2)]))
      Zcurr <- ZOSull1(xCurr, intKnots = intKnots, range.x = range.xVal)
      range.x.list[[jSpl]] <- attr(Zcurr, "range.x")
      intKnots.list[[jSpl]] <- attr(Zcurr, "intKnots")
      ncZspl_p <- c(ncZspl_p, ncol(Zcurr))
      Zspl_p <- cbind(Zspl_p, Zcurr)
      
    }
  }
  return(list(Z=Zspl_p,Knot_position=intKnots.list))
}

## Adding a O'Sullivan spline class and methods
# smooth.construct.os.smooth.spec<-function(object,data,knots) {
#   
#   x <- data[[object$term]]
#   num.knots <- min(25,round(length(x)/4)) ## number of knots
#   z.spline <- cal_z_kp(num.knots,as.matrix(x))
#   object$X<- cbind(1,x,z.spline$Z)
#   object$knots<- unlist(z.spline$Knot_position)
#   object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
#   object$S[[1]] <- OsplinePen(Boundary.knots=range(x), knots=unlist(z.spline$Knot_position), ord=2)
#   object$rank<- num.knots
#   object$bs.dim <- num.knots
#   class(object)<-"os.smooth"  # Give object a class
#   object
# }
# 
# Predict.matrix.os.smooth<-function(object,data) {
#   ## prediction method function for the `os' smooth class
#   x <- data[[object$term]]
#   num.knots <- min(25,round(length(x)/4)) ## number of knots
#   z.spline <- cal_z_kp(num.knots,as.matrix(x))
#   X <- cbind(1,x,z.spline$Z) # return the prediction matrix
#   X
# }