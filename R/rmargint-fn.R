#' Derivative of Tukey's bi-square loss function.
#'
#' This function evaluates the first derivative of Tukey's bi-square loss function.
#'
#' This function evaluates the first derivative of Tukey's bi-square loss function.
#'
#' @param r A vector of real numbers
#' @param k A positive tuning constant.
#'
#' @return A vector of the same length as \code{r}.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#'
#' @examples
#' x <- seq(-2, 2, length=10)
#' psi.tukey(r=x, k = 1.5)
#'
#' @export
#' @import stats graphics
#' @useDynLib rmargint, .registration = TRUE
psi.tukey <- function(r, k=4.685){
  u <- abs(r/k)
  w <- r*((1-u)*(1+u))^2
  w[u>1] <- 0
  return(w) 
}



# Tukey's weight function "Psi(r)/r"
psi.w <- function(r, k= 4.685){
  u <- abs(r/k)
  w <- ((1 + u) * (1 - u))^2
  w[u > 1] <- 0
  return(w)
}


#' Derivative of Huber's loss function.
#' 
#' This function evaluates the first derivative of Huber's loss function.
#'
#' This function evaluates the first derivative of Huber's loss function.
#' 
#' @param r A vector of real numbers.
#' @param k A positive tuning constant.
#' 
#' @return A vector of the same length as \code{r}.
#' 
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#' 
#' @examples
#' x <- seq(-2, 2, length=10)
#' psi.huber(r=x, k = 1.5)
#' 
#' @export
psi.huber <- function(r, k=1.345)
  pmin(k, pmax(-k, r))


#Huber's weight function "Psi(r)/r"
psi.huber.w <- function(r, k=1.345)
  pmin(1, k/abs(r))


#' Euclidean norm of a vector
#' 
#' This function calculates the Euclidean norm of a vector.
#' 
#' @param x A real vector.
#' 
#' @return The Euclidean norm of the input vector.
#' 
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#' 
#' @examples 
#' x <- seq(-2, 2, length=10)
#' my.norm.2(x)
#' 
#' @export
my.norm.2 <- function(x) sqrt(sum(x^2))


#' Epanechnikov kernel
#'
#' This function evaluates an Epanechnikov kernel
#'
#' This function evaluates an Epanechnikov kernel.
#'
#' @param x a vector of real numbers
#'
#' @return A vector of the same length as \code{x} where each entry is
#' \code{0.75 * (1 - x^2)} if \code{x < 1} and 0 otherwise.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#'
#' @examples
#' x <- seq(-2, 2, length=10)
#' k.epan(x)
#'
#' @export
k.epan<-function(x) {
  a <- 0.75*(1-x^2)
  tmp <- a*(abs(x)<=1)
  return(tmp)
}


#Order 2 kernel = Epanechnikov kernel
kernel2<-function(t){
  nucleo<- (3/4)*(1-t^2)*(abs(t)<=1)
  nucleo
}


#- Higher order kernels -#

#' Order 4 kernel
#' 
#' This function evaluates a kernel of order 4.
#' 
#' This function evaluates a kernel of order 4. A kernel L is a kernel of order 4 if it integrates 1, the integrals of u^j L(u) are 0 for 1 <= j < 4 (j integer) and the integral of u^4 L(u) is different from 0.
#' 
#' @param x A vector of real numbers.
#' 
#' @return A vector of the same length as \code{x} where each entry is \code{( 15/32 ) * ( 1 - x^2 ) * ( 3 - 7 * x^2 )} if \code{abs(x) < 1} and 0 otherwise.
#' 
#' @author Alejandra Martinez, \email{ale_m_martinez@hotmail.com}, Matias Salibian-Barrera
#' 
#' @examples 
#' x <- seq(-2,2,length=10)
#' kernel4(x)
#' 
#' @export
kernel4<-function(x) {
  a <- (15/32)*(1-x^2)*(3-7*x^2)
  tmp <- a*(abs(x)<= 1)
  return(tmp)
}


#' Order 6 kernel
#' 
#' This function evaluates a kernel of order 6.
#' 
#' This function evaluates a kernel of order 6. A kernel L is a kernel of order 6 if it integrates 1, the integrals of u^j L(u) are 0 for 1 <= j < 6 (j integer) and the integral of u^6 L(u) is different from 0.
#' 
#' 
#' @param x A vector of real numbers.
#' 
#' @return A vector of the same length as \code{x} where each entry is \code{( 105/256 ) * ( 1 - x^2 ) * ( 5 - 30 * x^2 + 33 * x^4 )} if \code{abs(x) < 1} and 0 otherwise.
#' 
#' @author Alejandra Martinez, \email{ale_m_martinez@hotmail.com}, Matias Salibian-Barrera
#' 
#' @examples 
#' x <- seq(-2,2,length=10)
#' kernel6(x)
#' 
#' @export
kernel6<-function(x) {
  a <- (105/256)*(1-x^2)*(5-30*x^2+33*x^4)
  tmp <- a*(abs(x)<=1)
  return(tmp)
}

#' Order 8 kernel
#' 
#' This function evaluates a kernel of order 8.
#' 
#' This function evaluates a kernel of order 8. A kernel L is a kernel of order 8 if it integrates 1, the integrals of u^j L(u) are 0 for 1 <= j < 8 (j integer) and the integral of u^8 L(u) is different from 0.
#' 
#' @param x A vector of real numbers.
#' 
#' @return A vector of the same length as \code{x} where each entry is \code{( 315/4096 ) * ( 1 - x^2 ) * ( 35 - 385 * x^2 + 1001 * x^4 - 715 * x^6 )} and 0 otherwise.
#' 
#' @author Alejandra Martinez, \email{ale_m_martinez@hotmail.com}, Matias Salibian-Barrera
#' 
#' @examples 
#' x <- seq(-2,2,length=10)
#' kernel8(x)
#' @export
kernel8<-function(x) {
  a <- (315/4096)*(1-x^2)*(35-385*x^2+1001*x^4-715*x^6)
  tmp <- a*(abs(x)<=1)
  return(tmp)
}

#' Order 10 kernel
#' 
#' This function evaluates a kernel of order 10. A kernel of order 10.
#' 
#' This function evaluates a kernel of order 10. A kernel L is a kernel of order 10 if it integrates 1, the integrals of u^j L(u) are 0 for 1 <= j < 10 (j integer) and the integral of u^10 L(u) is different from 0.
#' 
#' @param x A vector of real numbers.
#' 
#' @return A vector of the same length as \code{x} where each entry is \code{0.75 * ( 1 - x^2 ) * ( 315/128 - 105/32 * x^2 + 63/64 * x^4 - 3/32 * x^6 - 1/384 * x^8 )} and 0 otherwise.
#' 
#' @author Alejandra Martinez, \email{ale_m_martinez@hotmail.com}, Matias Salibian-Barrera
#' 
#' @examples 
#' x <- seq(-2,2,length=10)
#' kernel10(x)
#' 
#' @export
kernel10<-function(x) {
  a <- 0.75*(1-x^2)*(315/128-105/32*x^2+63/64*x^4-3/32*x^6-1/384*x^8)
  tmp <- a*(abs(x)<=1)
  return(tmp)
}

#' Classic marginal integration procedures for additive models
#' 
#' This function computes the standard marginal integration procedures for additive models.
#' 
#' This function computes three types of classical marginal integration procedures for additive models, that is, considering a squared loss function.
#' 
#' @param Xp Matrix (n by p) of explanatory variables.
#' @param yp  Vector of responses (missing values are allowed).
#' @param point Matrix of points where predictions will be computed and returned.
#' @param windows Vector or a squared matrix of bandwidths for the smoothing estimation procedure.
#' @param epsilon Convergence criterion.
#' @param prob Vector of robabilities of observing each response (n). Defaults to \code{NULL}.
#' @param type Three different type of estimators can be selected: type \code{'0'} (local constant on all the covariates), type \code{'1'} (local linear smoother on all the covariates), type \code{'alpha'} (local polynomial smoother only on the direction of interest).
#' @param degree Degree of the local polynomial smoother in the direction of interest when using the estimator of type \code{'alpha'}. Defaults to \code{NULL} for the case when using estimators of type \code{'0'} or \code{'1'}.
#' @param orderkernel Order of the kernel used in the nuisance directions when using the estimator of type \code{'alpha'}. Defaults to \code{2}.
#' @param qderivate If TRUE, it calculates \code{g^(q+1)/(q+1)!} for each component only for the type \code{'alpha'} method. Defaults to \code{FALSE}.
#' @param Qmeasure A matrix of points where the integration procedure ocurrs. Defaults to \code{NULL} for calcuting the integrals over the sample.
#' 
#' @return A list with the following components:
#' \item{mu}{Estimate for the intercept.}
#' \item{g.matrix}{Matrix of estimated additive components (n by p).}
#' \item{prediction }{Matrix of estimated additive components for the points listed in the argument point.}
#' \item{mul}{A vector of size p showing in each component the estimated intercept that considers only that direction of interest when using the type \code{'alpha'} method.}
#' \item{g.derivative }{Matrix of estimated derivatives of the additive components (only when qderivate is \code{TRUE}) (n by p).}
#' \item{prediction.derivate }{Matrix of estimated derivatives of the additive components for the points listed in the argument point (only when qderivate is \code{TRUE}).}
#' \item{Xp}{Matrix of explanatory variables.}
#' \item{yp}{Vector of responses.}
#' 
#' @references 
#' Chen R., Hardle W., Linton O.B. and Severance-Lossin E. (1996). Nonparametric estimation of additive separable regression models. Physica-Verlag HD, Switzerland.
#' Linton O. and Nielsen J. (1995). A kernel method of estimating structured nonparametric regression based on marginal integration. Biometrika, 82(1), 93-101.
#' Severance-Lossin E. and Sperlich S. (1999). Estimation of derivatives for additive separable models. Statistics, 33(3), 241-265.
#' Tjostheim D. and Auestad B. (1994). Nonparametric identification of nonlinear time series: Selecting significant lags. Journal of the American Statistical Association, 89(428), 1410-1430.
#' 
#' @author Alejandra Martinez, \email{ale_m_martinez@hotmail.com}, Matias Salibian-Barrera
#' 
#' @examples 
#' function.g1 <- function(x1) 24*(x1-1/2)^2-2
#' function.g2 <- function(x2) 2*pi*sin(pi*x2)-4
#' n <- 150
#' x1 <- runif(n)
#' x2 <- runif(n)
#' X <- cbind(x1, x2)
#' eps <- rnorm(n,0,sd=0.15)
#' regresion <- function.g1(x1) + function.g2(x2)
#' y <- regresion + eps
#' bandw <- matrix(0.25,2,2)
#' set.seed(8090)
#' nQ <- 80 
#' Qmeasure <- matrix(runif(nQ*2), nQ, 2)
#' fit.cl <- margint.cl(Xp=X, yp=y, windows=bandw, type='alpha', degree=1, Qmeasure=Qmeasure)
#' 
#' @export
margint.cl <- function(Xp, yp, point=NULL, windows, epsilon=1e-6, prob=NULL,
                       type='0', degree=NULL, qderivate=FALSE, orderkernel=2,
                       Qmeasure=NULL) {

  if(!is.null(dim(Xp))){
    if(type=='alpha'){
      if(is.null(degree)){
        return("Degree of local polynomial missing") #stop("Degree of local polynomial missing")
      }else{
        if( length(windows)==1 ){
          return("Windows should be a vector o a matrix") #stop("Windows should be a vector o a matrix")
        }
      }
    }else{
      if( (is.matrix(windows)) | (length(windows)==1)  ){
        return("Windows should be a vector") #stop("Windows should be a vector")
      }
    }
  }else{
    if(!is.null(dim(windows))){
      return("Windows should be a number") #stop("Windows should be a number")
    }
  }


  n <- length(yp)
  Xp <- as.matrix(Xp)
  q <- dim(Xp)[2]
  corte <- 10*epsilon
  punto <- point

  # Remove observations with missing responses
  yy <- yp
  XX <- Xp
  yp <- yp[ tmp<-!is.na(yp) ]
  Xp <- Xp[tmp, , drop=FALSE]
  n.miss <- length(yp)
  if(is.null(prob)){prob <- rep(1,n.miss)
  } else {
    prob <- prob[tmp]
  }
  if(dim(t(as.matrix(windows)))[2]!=q){return("Error Dimension of Bandwidths")}


  #Initializations
  puntoj <- rep(0,q)
  pesosi <- rep(1,n.miss)

  # If a Qmeasure is provided
  if(!is.null(Qmeasure)){

    grilla <- rbind(Qmeasure, XX)
    nQ <- dim(grilla)[1]

    alphal.aux <- matrix(0,nQ,q)
    g.matriz.bis <- matrix(0,nQ,q)
    g.matriz <- matrix(0,n,q)
    aux.b <- rep(0,nQ)
    alpha.aux <- rep(0, nQ)
    nq <- dim(Qmeasure)[1]
    aa <- rep(0,nq)

    if(qderivate){
      g.derivate.bis <- matrix(0,nQ,q)
      g.derivate <- matrix(0,n,q)
      aa.derivate <- rep(0,nq)
    }else{
      g.derivate <- NULL
    }


    for(k in 1:nQ){
      for(j in 1:q){
        for(i in 1:nq){
          puntoj <- Qmeasure[i,]
          puntoj[j]<- grilla[k,j]
          if(type=='0'){
            aa[i] <- .C("kernel_cl_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                        as.double(yp), as.double(windows), as.double(epsilon), as.double(prob),salida=as.double(0) )$salida
            if(i==k){alpha.aux[k] <- aa[i]}
          }
          if(type=='1'){
            aa[i] <- .C("kernel_cl_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                        as.double(yp), as.double(windows), as.double(epsilon), as.double(prob),salida=as.double(rep(0, q+1)))$salida[1]
            if(i==k){alpha.aux[k] <- aa[i]}
          }
          if(type=='alpha'){
            if(is.null(dim(windows))){
              windows <- t(matrix(windows,q,q))
            }
            AUX <- .C("kernel_cl_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                      as.double(yp), as.double(windows[j,]), as.double(epsilon), as.double(prob), salida=as.double(rep(0, degree+1)) )$salida
            aa[i] <- AUX[1]
            if(qderivate){
              aa.derivate[i] <- AUX[degree+1]
            }
            if(i==k){alphal.aux[k,j] <- aa[i]}
          }
        }
        g.matriz.bis[k,j]<-mean(aa,na.rm=TRUE)
        if(qderivate){
          g.derivate.bis[k,j]<-mean(aa.derivate,na.rm=TRUE)
        }
      }
    }

    alphal <- NULL
    if(type=='alpha'){
      alpha.aux <- colMeans(alphal.aux[1:nq,],na.rm=TRUE)
      alphal <- alpha.aux
    }

    alpha <- mean(alpha.aux,na.rm=TRUE)
    g.matriz <- g.matriz.bis[(nq+1):nQ,] - alpha
    if(qderivate){
      g.derivate <- g.derivate.bis[(nq+1):nQ,]
    }
  }else{ #That is, if a Qmeasure is not provided

    alphal.aux <- matrix(0,n,q)
    g.matriz <- matrix(0,n,q)
    if(qderivate){
      g.derivate <- matrix(0,n,q)
      A.derivate <- rep(0,n)
    }else{
      g.derivate <- NULL
    }
    aux.b <- rep(0,n)
    aa <- rep(0,n)
    alpha.aux <- rep(0, n)

    for(k in 1:n){
      for(j in 1:q){
        for(i in 1:n){
          puntoj <- XX[i,]
          puntoj[j]<-XX[k,j]
          if(type=='0'){
            aa[i] <- .C("kernel_cl_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                        as.double(yp), as.double(windows), as.double(epsilon), as.double(prob),salida=as.double(0) )$salida
            if(i==k){alpha.aux[k] <- aa[i]}
          }
          if(type=='1'){
            aa[i] <- .C("kernel_cl_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                        as.double(yp), as.double(windows), as.double(epsilon), as.double(prob),salida=as.double(rep(0, q+1)))$salida[1]
            if(i==k){alpha.aux[k] <- aa[i]}
          }
          if(type=='alpha'){
            if(is.null(dim(windows))){
              windows <- t(matrix(windows,q,q))
            }
            AUX <- .C("kernel_cl_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                      as.double(yp), as.double(windows[j,]), as.double(epsilon), as.double(prob), salida=as.double(rep(0, degree+1)) )$salida
            aa[i] <- AUX[1]
            if(qderivate){
              A.derivate[i] <- AUX[degree+1]
            }
            if(i==k){alphal.aux[k,j] <- aa[i]}
          }
        }
        g.matriz[k,j]<-mean(aa,na.rm=TRUE)
        if(qderivate){
          g.derivate[k,j]<-mean(A.derivate,na.rm=TRUE)
        }
      }
    }

    alphal <- NULL
    if(type=='alpha'){
      alpha.aux <- colMeans(alphal.aux,na.rm=TRUE)
      alphal <- alpha.aux
    }

    alpha <- mean(alpha.aux,na.rm=TRUE)
    g.matriz <- g.matriz - alpha
  }

  #Predictions:

  prediccion <- NULL

  #If a Qmeasure is provided
  if(!is.null(Qmeasure)){
    aa <- rep(0,nq)
    aa.deri <- rep(0,nq)
    if(!is.null(punto)) {
      if(is.null(dim(punto))) {
        prediccion <- mpunto <- t(as.matrix(punto))
        if(qderivate){
          prediccion.deri <- prediccion
        }
      } else {
        prediccion <- mpunto <- punto
        if(qderivate){
          prediccion.deri <- prediccion
        }
      }
      np <- dim(mpunto)[1]
      for(k in 1:np){
        for(j in 1:q){
          for(i in 1:nq){
            puntoj <- Qmeasure[i,]
            puntoj[j] <- mpunto[k,j]

            if(type=='0'){
              aa[i] <- .C("kernel_cl_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(windows), as.double(epsilon), as.double(prob),salida=as.double(0) )$salida
            }
            if(type=='1'){
              aa[i] <- .C("kernel_cl_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(windows), as.double(epsilon), as.double(prob),salida=as.double(rep(0, q+1)))$salida[1]
            }

            if(type=='alpha'){
              AUX <- .C("kernel_cl_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                        as.double(yp), as.double(windows[j,]), as.double(epsilon), as.double(prob), salida=as.double(rep(0, degree+1)) )$salida
              aa[i] <- AUX[1]
              if(qderivate){
                aa.deri[i] <- AUX[degree+1]
              }
            }
          }
          prediccion[k,j] <- mean(aa,na.rm=TRUE)-alpha
          if(qderivate){
            prediccion.deri[k,j] <- mean(aa.deri,na.rm=TRUE)
          }
        }
      }
    } #end if
  }else{ #If a Qmeasure is not provided
    aa <- rep(0,n)
    aa.deri <- rep(0,n)
    if(!is.null(punto)){
      if(is.null(dim(punto))){
        prediccion <- mpunto <- t(as.matrix(punto))
        if(qderivate){
          prediccion.deri <- prediccion
        }
      } else {
        prediccion <- mpunto <- punto
        if(qderivate){
          prediccion.deri <- prediccion
        }
      }
      np <- dim(mpunto)[1]
      for(k in 1:np){
        for(j in 1:q){
          for(i in 1:n){
            puntoj <- XX[i,]
            puntoj[j] <- mpunto[k,j]

            if(type=='0'){
              aa[i] <- .C("kernel_cl_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(windows), as.double(epsilon), as.double(prob),salida=as.double(0) )$salida
            }
            if(type=='1'){
              aa[i] <- .C("kernel_cl_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(windows), as.double(epsilon), as.double(prob),salida=as.double(rep(0, q+1)))$salida[1]
            }

            if(type=='alpha'){
              AUX <- .C("kernel_cl_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                        as.double(yp), as.double(windows[j,]), as.double(epsilon), as.double(prob), salida=as.double(rep(0, degree+1)) )$salida
              aa[i] <- AUX[1]
              if(qderivate){
                aa.deri[i] <- AUX[degree+1]
              }
            }
          }
          prediccion[k,j] <- mean(aa,na.rm=TRUE)-alpha
          if(qderivate){
            prediccion.deri[k,j] <- mean(aa.deri,na.rm=TRUE)
          }
        }
      }
    } #End of if
  } #End of else

  if(!is.null(point)){
    if(type=='alpha'){
      if(!qderivate){
        object <- list(mu=alpha,g.matrix=g.matriz, prediction=prediccion, mul=alphal, Xp=Xp, yp=yp)
        class(object) <- c("margint.cl", "margint", "list")
        return(object)
      } else {
        object <- list(mu=alpha,g.matrix=g.matriz, prediction=prediccion, mul=alphal,g.derivate=g.derivate, prediction.derivate=prediccion.deri, Xp=Xp, yp=yp)
        class(object) <- c("margint.cl", "margint", "list")
        return(object)
      }
    } else {
      object <- list(mu=alpha,g.matrix=g.matriz, prediction=prediccion, Xp=Xp, yp=yp)
      class(object) <- c("margint.cl", "margint", "list")
      return(object)
    }
  } else {
    if(type=='alpha'){
      if(!qderivate){
        object <- list(mu=alpha,g.matrix=g.matriz, mul=alphal, Xp=Xp, yp=yp)
        class(object) <- c("margint.cl", "margint", "list")
        return(object)
      } else {
        object <- list(mu=alpha,g.matrix=g.matriz, mul=alphal,g.derivate=g.derivate, Xp=Xp, yp=yp)
        class(object) <- c("margint.cl", "margint", "list")
        return(object)
      }
    } else {
      object <- list(mu=alpha,g.matrix=g.matriz, Xp=Xp, yp=yp)
      class(object) <- c("margint.cl", "margint", "list")
      return(object)
    }
  }
}

#' Robust marginal integration procedures for additive models
#' 
#' This function computes robust marginal integration procedures for additive models.
#' 
#' This function computes three types of robust marginal integration procedures for additive models.
#' 
#' @param Xp Matrix (n by p) of explanatory variables.
#' @param yp  Vector of responses (missing values are allowed).
#' @param point Matrix of points where predictions will be computed and returned.
#' @param windows Vector or a squared matrix of bandwidths for the smoothing estimation procedure.
#' @param prob Probabilities of observing each response (n). Defaults to \code{NULL}.
#' @param sigma.hat Estimate of the residual standard error. If \code{NULL} we use the mad of the residuals obtained with local medians.
#' @param win.sigma Vector of bandwidths for estimating sigma.hat. If \code{NULL} it uses the argument windows if it is a vector or its diagonal if it is a matrix.
#' @param epsilon Convergence criterion.
#' @param type Three different type of estimators can be selected: type \code{'0'} (local constant on all the covariates), type \code{'1'} (local linear smoother on all the covariates), type \code{'alpha'} (local polynomial smoother only on the direction of interest).
#' @param degree Degree of the local polynomial smoother in the direction of interest when using the estimator of type \code{'alpha'}. Defaults to \code{NULL} for the case when using estimators of type \code{'0'} or \code{'1'}.
#' @param typePhi One of either \code{'Tukey'} or \code{'Huber'}.
#' @param k.h Tuning constant for a Huber-type loss function. Defaults to \code{1.345}.
#' @param k.t Tuning constant for a Tukey-type loss function. Defaults to \code{4.685}.
#' @param max.it Maximum number of iterations for the algorithm.
#' @param qderivate If TRUE, it calculates \code{g^(q+1)/(q+1)!} for each component only for the type \code{'alpha'} method. Defaults to \code{FALSE}.
#' @param orderkernel Order of the kernel used in the nuisance directions when using the estimator of type \code{'alpha'}. Defaults to \code{2}.
#' @param Qmeasure A matrix of points where the integration procedure ocurrs. Defaults to \code{NULL} for calcuting the integrals over the sample.
#' 
#' @return A list with the following components:
#' \item{mu }{Estimate for the intercept.}
#' \item{g.matrix }{Matrix of estimated additive components (n by p).}
#' \item{sigma.hat }{Estimate of the residual standard error.}
#' \item{prediction }{Matrix of estimated additive components for the points listed in the argument point.}
#' \item{mul }{A vector of size p showing in each component the estimated intercept that considers only that direction of interest when using the type \code{'alpha'} method.}
#' \item{g.derivative }{Matrix of estimated derivatives of the additive components (only when qderivate is \code{TRUE}) (n by p).}
#' \item{prediction.derivate }{Matrix of estimated derivatives of the additive components for the points listed in the argument point (only when qderivate is \code{TRUE}).}
#' \item{Xp}{Matrix of explanatory variables.}
#' \item{yp}{Vector of responses.}
#'
#' @author Alejandra Martinez, \email{ale_m_martinez@hotmail.com}, Matias Salibian-Barrera
#' 
#' @references Boente G. and Martinez A. (2017). Marginal integration M-estimators for additive models. TEST, 26(2), 231-260. https://doi.org/10.1007/s11749-016-0508-0
#' 
#' @examples 
#' function.g1 <- function(x1) 24*(x1-1/2)^2-2
#' function.g2 <- function(x2) 2*pi*sin(pi*x2)-4
#' set.seed(140)
#' n <- 150
#' x1 <- runif(n)
#' x2 <- runif(n)
#' X <- cbind(x1, x2)
#' eps <- rnorm(n,0,sd=0.15)
#' regresion <- function.g1(x1) + function.g2(x2)
#' y <- regresion + eps
#' bandw <- matrix(0.25,2,2)
#' set.seed(8090)
#' nQ <- 80 
#' Qmeasure <- matrix(runif(nQ*2), nQ, 2)
#' fit.rob <- margint.rob(Xp=X, yp=y, windows=bandw, type='alpha', degree=1, Qmeasure=Qmeasure) 
#' 
#' @export
margint.rob <- function(Xp, yp, point=NULL, windows, prob=NULL, sigma.hat=NULL,
                        win.sigma=NULL, epsilon=1e-06, type='0', degree=NULL, typePhi='Huber',
                        k.h=1.345, k.t = 4.685, max.it=20, qderivate=FALSE, orderkernel=2,
                        Qmeasure=NULL){

  
  if(!is.null(dim(Xp))){
    if(type=='alpha'){
      if(is.null(degree)){
        return("Degree of local polynomial missing") #stop("Degree of local polynomial missing")
      }else{
        if( length(windows)==1 ){
          return("Windows should be a vector or a matrix") #stop("Windows should be a vector or a matrix")
        }
      }
    }else{
      if( (is.matrix(windows)) | (length(windows)==1) ){
        return("Windows should be a vector") #stop("Windows should be a vector")
      }
    }
  }else{
    if(!is.null(dim(windows))){
      return("Windows should be a number") #stop("Windows should be a number")
    }
  }

  n <- length(yp)
  Xp <- as.matrix(Xp)
  q <- dim(Xp)[2]
  corte <- 10*epsilon
  punto <- point

  if(is.null(win.sigma)){
    if(is.null(dim(windows))){
      win.sigma <- windows
    }else{
      win.sigma <- diag(windows)
    }
  }

  # Remove observations with missing responses
  yy<-yp
  XX <- Xp
  yp <- yp[ tmp<-!is.na(yp) ]
  Xp <- Xp[tmp, ,drop=FALSE]
  n.miss <- length(yp)
  if(is.null(prob)){prob <- rep(1,n.miss)
  }else{
    prob <- prob[tmp]
  }
  if(dim(t(as.matrix(windows)))[2]!=q){return("Error Dimension of Bandwidths")}

  #Initializations
  puntoj <- rep(0,q)
  aa <- rep(0,n)
  pesosi <- rep(1,n.miss)

  # Estimate residual standard error
  if(is.null(sigma.hat)){
    ab <- rep(0,n.miss)
    for(i in 1:n.miss){
      xtildebis <- scale(Xp, center=Xp[i,], scale=win.sigma)
      a <- matrix(as.numeric(abs(xtildebis) < 1), n.miss, q)
      a <- apply(a, 1, prod)
      a[ a == 0 ] <- NA
      ab[i] <- median( a*yp, na.rm = TRUE)
    }
    sigma.hat <- mad(yp - ab)
    if( sigma.hat < 1e-10 ){sigma.hat <- 1e-10} # sigma.hat <- sd(yp-ab,na.rm=TRUE)
  }


  #If a Qmeasure is provided
  if(!is.null(Qmeasure)){

    grilla <- rbind(Qmeasure, XX)
    nQ <- dim(grilla)[1]

    alphal.aux <- matrix(0,nQ,q)
    g.matriz.bis <- matrix(0,nQ,q)
    g.matriz <- matrix(0,n,q)
    aux.b <- rep(0,nQ)
    alpha.aux <- rep(0, nQ)
    nq <- dim(Qmeasure)[1]
    aa <- rep(0,nq)
    if(qderivate){
      g.derivate.bis <- matrix(0,nQ,q)
      g.derivate <- matrix(0,n,q)
      aa.derivate <- rep(0,nq)
    }else{
      g.derivate <- NULL
    }
    for(k in 1:nQ){
      for(j in 1:q){
        for(i in 1:nq){
          puntoj <- Qmeasure[i,]
          puntoj[j]<-grilla[k,j]

          #Inicializo el mu
          if(type=='0'){
            mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                         as.double(yp), as.double(windows),
                         as.double(prob), salida=as.double(0) )$salida
            if(typePhi=='Huber'){
              aa[i] <- .C("kernel_huber_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(mu.ini), as.double(windows),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(0) )$salida
            }
            if(typePhi=='Tukey'){
              aa[i] <- .C("kernel_tukey_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(mu.ini), as.double(windows),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(0) )$salida
            }
            if(i==k){alpha.aux[k] <- aa[i]}
          }
          if(type=='1'){
            mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                         as.double(yp), as.double(windows),
                         as.double(prob), salida=as.double(0) )$salida
            beta.ini <- rep(0,q+1)
            beta.ini[1] <- mu.ini
            if(typePhi=='Huber'){
              aa[i] <- .C("kernel_huber_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(beta.ini), as.double(windows),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(rep(0, q+1)) )$salida[1]

            }
            if(typePhi=='Tukey'){
              aa[i] <- .C("kernel_tukey_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(beta.ini), as.double(windows),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(rep(0, q+1)) )$salida[1]
            }
            if(i==k){alpha.aux[k] <- aa[i]}
          }
          if(type=='alpha'){
            mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                         as.double(yp), as.double(windows),
                         as.double(prob), salida=as.double(0) )$salida
            beta.ini <- rep(0,degree+1)
            beta.ini[1] <- mu.ini
            if(is.null(dim(windows))){
              windows <- t(matrix(windows,q,q))
            }
            if(typePhi=='Huber'){
              AUX <- .C("kernel_huber_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                        as.double(yp), as.double(beta.ini), as.double(windows[j,]),
                        as.double(epsilon), as.double(sigma.hat),
                        as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(rep(0, degree+1)) )$salida
              aa[i] <- AUX[1]
              if(qderivate){
                aa.derivate[i] <- AUX[degree+1]
              }
            }
            if(typePhi=='Tukey'){
              AUX <- .C("kernel_tukey_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                        as.double(yp), as.double(beta.ini), as.double(windows[j,]),
                        as.double(epsilon), as.double(sigma.hat),
                        as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(rep(0, degree+1)) )$salida
              aa[i] <- AUX[1]
              if(qderivate){
                aa.derivate[i] <- AUX[degree+1]
              }
            }
            if(i==k){
              alphal.aux[i,j] <- aa[i]
            }
          }
        }
        g.matriz.bis[k,j]<-mean(aa,na.rm=TRUE)
        if(qderivate){
          g.derivate.bis[k,j]<-mean(aa.derivate,na.rm=TRUE)
        }
      }
    }

    alphal <- NULL
    if(type=='alpha'){
      alpha.aux <- colMeans(alphal.aux[1:nq,],na.rm=TRUE)
      alphal <- alpha.aux
    }
    alpha <- mean(alpha.aux,na.rm=TRUE)
    g.matriz <- g.matriz.bis[(nq+1):nQ,] - alpha
    if(qderivate){
      g.derivate <- g.derivate.bis[(nq+1):nQ,]
    }
  }else{ #If no Qmeasure is provided
    alpha.aux <- rep(0, n)
    g.matriz <- matrix(0,n,q)
    if(qderivate){
      g.derivate <- matrix(0,n,q)
      aa.derivate <- rep(0,n)
    }else{
      g.derivate <- NULL
    }
    alphal.aux <- matrix(0,n,q)

    for(k in 1:n){
      for(j in 1:q){
        for(i in 1:n){
          puntoj <- XX[i,]
          puntoj[j]<-XX[k,j]

          #Inicializo el mu
          if(type=='0'){
            mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                         as.double(yp), as.double(windows),
                         as.double(prob), salida=as.double(0) )$salida
            if(typePhi=='Huber'){
              aa[i] <- .C("kernel_huber_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(mu.ini), as.double(windows),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(0) )$salida
            }
            if(typePhi=='Tukey'){
              aa[i] <- .C("kernel_tukey_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(mu.ini), as.double(windows),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(0) )$salida
            }
            if(i==k){alpha.aux[k] <- aa[i]}
          }
          if(type=='1'){
            mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                         as.double(yp), as.double(windows),
                         as.double(prob), salida=as.double(0) )$salida
            beta.ini <- rep(0,q+1)
            beta.ini[1] <- mu.ini
            if(typePhi=='Huber'){
              aa[i] <- .C("kernel_huber_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(beta.ini), as.double(windows),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(rep(0, q+1)) )$salida[1]

            }
            if(typePhi=='Tukey'){
              aa[i] <- .C("kernel_tukey_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                          as.double(yp), as.double(beta.ini), as.double(windows),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(rep(0, q+1)) )$salida[1]
            }
            if(i==k){alpha.aux[k] <- aa[i]}
          }
          if(type=='alpha'){
            mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                         as.double(yp), as.double(windows),
                         as.double(prob), salida=as.double(0) )$salida
            beta.ini <- rep(0,degree+1)
            beta.ini[1] <- mu.ini
            if(is.null(dim(windows))){
              windows <- t(matrix(windows,q,q))
            }
            if(typePhi=='Huber'){
              AUX <- .C("kernel_huber_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                        as.double(yp), as.double(beta.ini), as.double(windows[j,]),
                        as.double(epsilon), as.double(sigma.hat),
                        as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(rep(0, degree+1)) )$salida
              aa[i] <- AUX[1]
              if(qderivate){
                aa.derivate[i] <- AUX[degree+1]
              }
            }
            if(typePhi=='Tukey'){
              AUX <- .C("kernel_tukey_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                        as.double(yp), as.double(beta.ini), as.double(windows[j,]),
                        as.double(epsilon), as.double(sigma.hat),
                        as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(rep(0, degree+1)) )$salida
              aa[i] <- AUX[1]
              if(qderivate){
                aa.derivate[i] <- AUX[degree+1]
              }
            }
            if(i==k){
              alphal.aux[i,j] <- aa[i]
            }
          }
        }
        g.matriz[k,j]<-mean(aa,na.rm=TRUE)
        if(qderivate){
          g.derivate[k,j]<-mean(aa.derivate,na.rm=TRUE)
        }
      }
    }
    alphal <- NULL
    if(type=='alpha'){
      alpha.aux <- colMeans(alphal.aux,na.rm=TRUE)
      alphal <- alpha.aux
    }
    alpha <- mean(alpha.aux,na.rm=TRUE)
    g.matriz <- g.matriz - alpha
  }#End of else

  #Predictions:

  prediccion <- NULL

  #If a Qmeasure is provided
  if(!is.null(Qmeasure)){
    aa <- rep(0,nq)
    aa.deri <- rep(0,nq)
    if(!is.null(punto)) {
      if(is.null(dim(punto))){
        prediccion <- mpunto <- t(as.matrix(punto))
        if(qderivate){
          prediccion.deri <- prediccion
        }
      } else {
        prediccion <- mpunto <- punto
        if(qderivate){
          prediccion.deri <- prediccion
        }
      }
      np <- dim(mpunto)[1]
      for(k in 1:np){
        for(j in 1:q){
          for(i in 1:nq){
            puntoj <- grilla[i,]
            puntoj[j] <- mpunto[k,j]
            if(type=='0'){
              #Inicializo el mu
              mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                           as.double(yp), as.double(windows),
                           as.double(prob), salida=as.double(0) )$salida

              if(typePhi=='Huber'){
                aa[i] <- .C("kernel_huber_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                            as.double(yp), as.double(mu.ini), as.double(windows),
                            as.double(epsilon), as.double(sigma.hat),
                            as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(0) )$salida
              }
              if(typePhi=='Tukey'){
                aa[i] <- .C("kernel_tukey_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                            as.double(yp), as.double(mu.ini), as.double(windows),
                            as.double(epsilon), as.double(sigma.hat),
                            as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(0) )$salida
              }
            }
            if(type=='1'){
              mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                           as.double(yp), as.double(windows),
                           as.double(prob), salida=as.double(0) )$salida
              beta.ini <- rep(0,q+1)
              beta.ini[1] <- mu.ini
              if(typePhi=='Huber'){
                aa[i] <- .C("kernel_huber_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                            as.double(yp), as.double(beta.ini), as.double(windows),
                            as.double(epsilon), as.double(sigma.hat),
                            as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(rep(0, q+1)) )$salida[1]
              }
              if(typePhi=='Tukey'){
                aa[i] <- .C("kernel_tukey_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                            as.double(yp), as.double(beta.ini), as.double(windows),
                            as.double(epsilon), as.double(sigma.hat),
                            as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(rep(0, q+1)) )$salida[1]
              }
            }
            if(type=='alpha'){
              mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                           as.double(yp), as.double(windows),
                           as.double(prob), salida=as.double(0) )$salida
              beta.ini <- rep(0,degree+1)
              beta.ini[1] <- mu.ini
              if(typePhi=='Huber'){
                AUX <- .C("kernel_huber_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                          as.double(yp), as.double(beta.ini), as.double(windows[j,]),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(rep(0, degree+1)) )$salida
                aa[i] <- AUX[1]
                if(qderivate){
                  aa.deri[i] <- AUX[degree+1]
                }
              }
              if(typePhi=='Tukey'){
                AUX <- .C("kernel_tukey_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                          as.double(yp), as.double(beta.ini), as.double(windows[j,]),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(rep(0, degree+1)) )$salida
                aa[i] <- AUX[1]
                if(qderivate){
                  aa.deri[i] <- AUX[degree+1]
                }
              }
            }
          }
          prediccion[k,j] <- mean(aa,na.rm=TRUE)-alpha
          if(qderivate){
            prediccion.deri[k,j] <- mean(aa.deri,na.rm=TRUE)
          }
        }
      }
    }
  }else{ #If no Qmeasure is provided
    aa <- rep(0,n)
    aa.deri <- rep(0,n)
    if(!is.null(punto)) {
      if(is.null(dim(punto))){
        prediccion <- mpunto <- t(as.matrix(punto))
        if(qderivate){
          prediccion.deri <- prediccion
        }
      } else {
        prediccion <- mpunto <- punto
        if(qderivate){
          prediccion.deri <- prediccion
        }
      }
      np <- dim(mpunto)[1]
      for(k in 1:np){
        for(j in 1:q){
          for(i in 1:n){
            puntoj <- XX[i,]
            puntoj[j] <- mpunto[k,j]
            if(type=='0'){
              #Inicializo el mu
              mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                           as.double(yp), as.double(windows),
                           as.double(prob), salida=as.double(0) )$salida

              if(typePhi=='Huber'){
                aa[i] <- .C("kernel_huber_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                            as.double(yp), as.double(mu.ini), as.double(windows),
                            as.double(epsilon), as.double(sigma.hat),
                            as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(0) )$salida
              }
              if(typePhi=='Tukey'){
                aa[i] <- .C("kernel_tukey_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                            as.double(yp), as.double(mu.ini), as.double(windows),
                            as.double(epsilon), as.double(sigma.hat),
                            as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(0) )$salida
              }
            }
            if(type=='1'){
              mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                           as.double(yp), as.double(windows),
                           as.double(prob), salida=as.double(0) )$salida
              beta.ini <- rep(0,q+1)
              beta.ini[1] <- mu.ini
              if(typePhi=='Huber'){
                aa[i] <- .C("kernel_huber_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                            as.double(yp), as.double(beta.ini), as.double(windows),
                            as.double(epsilon), as.double(sigma.hat),
                            as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(rep(0, q+1)) )$salida[1]
              }
              if(typePhi=='Tukey'){
                aa[i] <- .C("kernel_tukey_lin_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                            as.double(yp), as.double(beta.ini), as.double(windows),
                            as.double(epsilon), as.double(sigma.hat),
                            as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(rep(0, q+1)) )$salida[1]
              }
            }
            if(type=='alpha'){
              mu.ini <- .C("ini_mu_pos_multi", puntoj, as.double(Xp), as.integer(q), as.integer(n.miss),
                           as.double(yp), as.double(windows),
                           as.double(prob), salida=as.double(0) )$salida
              beta.ini <- rep(0,degree+1)
              beta.ini[1] <- mu.ini
              if(typePhi=='Huber'){
                AUX <- .C("kernel_huber_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                          as.double(yp), as.double(beta.ini), as.double(windows[j,]),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(rep(0, degree+1)) )$salida
                aa[i] <- AUX[1]
                if(qderivate){
                  aa.deri[i] <- AUX[degree+1]
                }
              }
              if(typePhi=='Tukey'){
                AUX <- .C("kernel_tukey_alpha_multi", puntoj, as.double(Xp), as.integer(j), as.integer(degree), as.integer(q), as.integer(n.miss), as.integer(orderkernel),
                          as.double(yp), as.double(beta.ini), as.double(windows[j,]),
                          as.double(epsilon), as.double(sigma.hat),
                          as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(rep(0, degree+1)) )$salida
                aa[i] <- AUX[1]
                if(qderivate){
                  aa.deri[i] <- AUX[degree+1]
                }
              }
            }
          }
          prediccion[k,j] <- mean(aa,na.rm=TRUE)-alpha
          if(qderivate){
            prediccion.deri[k,j] <- mean(aa.deri,na.rm=TRUE)
          }
        }
      }
    }
  }


  if(!is.null(point)){
    if(type=='alpha'){
      if(!qderivate){
        object <- list(mu=alpha, g.matrix=g.matriz, sigma.hat=sigma.hat, prediction=prediccion, mul=alphal, Xp=Xp, yp=yp)
        class(object) <- c("margint.rob", "margint", "list")
        return(object)
      } else {
        object <- list(mu=alpha, g.matrix=g.matriz, sigma.hat=sigma.hat, prediction=prediccion, mul=alphal, g.derivate=g.derivate, prediction.derivate=prediccion.deri, Xp=Xp, yp=yp)
        class(object) <- c("margint.rob", "margint", "list")
        return(object)
      }
    } else {
      object <- list(mu=alpha, g.matrix=g.matriz, sigma.hat=sigma.hat, prediction=prediccion, Xp=Xp, yp=yp)
      class(object) <- c("margint.rob", "margint", "list")
      return(object)
    }
  } else {
    if(type=='alpha'){
      if(!qderivate){
        object <- list(mu=alpha, g.matrix=g.matriz, sigma.hat=sigma.hat, mul=alphal, Xp=Xp, yp=yp)
        class(object) <- c("margint.rob", "margint", "list")
        return(object)
      } else {
        object <- list(mu=alpha, g.matrix=g.matriz, sigma.hat=sigma.hat, mul=alphal,g.derivate=g.derivate, Xp=Xp, yp=yp)
        class(object) <- c("margint.rob", "margint", "list")
        return(object)
      }
    } else {
      object <- list(mu=alpha, g.matrix=g.matriz, sigma.hat=sigma.hat, Xp=Xp, yp=yp)
      class(object) <- c("margint.rob", "margint", "list")
      return(object)
    }
  }
}


#- S3 Methods -#

#' Residuals for objects of class \code{margint}
#'
#' This function returns the residuals of the fitted additive model using one of the three
#' classical or robust marginal integration estimators, as computed with \code{\link{margint.cl}} or
#' \code{\link{margint.rob}}.
#'
#' @param object an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.
#' @param ... additional other arguments. Currently ignored.
#'
#' @return A vector of residuals.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
residuals.margint <- function(object, ...){
  return( object$yp - rowSums(object$g.matrix) -object$mu )
}


#' Fitted values for objects of class \code{margint}
#'
#' This function returns the fitted values given the covariates of the original sample under an additive model using a classical or robust marginal integration procedure estimator computed with \code{\link{margint.cl}} or \code{\link{margint.rob}}.
#'
#' @param object an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.
#' @param ... additional other arguments. Currently ignored.
#'
#' @return A vector of fitted values.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
predict.margint <- function(object, ...){
  return( rowSums(object$g.matrix) + object$mu )
}

#' Diagnostic plots for objects of class \code{margint}
#'
#' Plot method for class \code{margint}.
#'
#' @param x an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.
#' @param derivative if TRUE, it plots the q-th derivatives. Defaults to FALSE.
#' @param which vector of indices of explanatory variables for which partial residuals plots will be generated. Defaults to all available explanatory variables.
#' @param ask logical value. If \code{TRUE}, the graphical device will prompt before going to the next page/screen of output.
#' @param ... additional other arguments.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#' 
#' @examples 
#' function.g1 <- function(x1) 24*(x1-1/2)^2-2
#' function.g2 <- function(x2) 2*pi*sin(pi*x2)-4
#' set.seed(140)
#' n <- 150
#' x1 <- runif(n)
#' x2 <- runif(n)
#' X <- cbind(x1, x2)
#' eps <- rnorm(n,0,sd=0.15)
#' regresion <- function.g1(x1) + function.g2(x2)
#' y <- regresion + eps
#' bandw <- matrix(0.25,2,2)
#' set.seed(8090)
#' nQ <- 80 
#' Qmeasure <- matrix(runif(nQ*2), nQ, 2)
#' fit.rob <- margint.rob(Xp=X, yp=y, windows=bandw, type='alpha', degree=1, Qmeasure=Qmeasure)
#' plot(fit.rob, which=1)
#' 
#' @export
plot.margint <- function(x, derivative=FALSE, which=1:np, ask=FALSE,...){
  object <- x
  Xp <- object$Xp
  np <- dim(Xp)[2]
  opar <- par(ask=ask)
  on.exit(par(opar))
  these <- rep(FALSE, np)
  these[ which ] <- TRUE
  if(!derivative){
    for(i in 1:np) {
      if(these[i]) {
        ord <- order(Xp[,i])
        if (is.null(colnames(Xp)) ){
          x_name <- bquote(paste('x')[.(i)])
          #paste("x",i,sep="")
        } else {
          x_name <- colnames(Xp)[i]
        }
        y_name <- bquote(paste(hat('g')[.(i)]))
        res <- object$yp - rowSums(object$g.matrix[,-i, drop=FALSE])-object$mu
        lim_cl <- c(min(res), max(res))
        plot(Xp[,i], res, pch=20,col='gray45',main="",xlab=x_name,ylab=y_name, ylim=lim_cl,cex.lab=0.8)
        lines(Xp[ord,i],object$g.matrix[ord,i],lwd=3)
      }
    }
  }else{
    for(i in 1:np) {
      if(these[i]) {
        ord <- order(Xp[,i])
        if (is.null(colnames(Xp)) ){
          x_name <- bquote(paste('x')[.(i)])
          #paste("x",i,sep="")
        } else {
          x_name <- colnames(Xp)[i]
        }
        y_name <- bquote(paste('g'[.(i)]^('q')))
        plot(Xp[ord,i], object$g.derivate[ord,i], type='l',lwd=3, main="",xlab=x_name,ylab=y_name,cex.lab=0.8)
      }
    }
  }
}

#' Summary for additive models fits using a marginal integration procedure
#'
#' Summary method for class \code{margint}.
#'
#' This function returns the estimation of the intercept and also the five-number summary and the mean of the residuals for both classical and robust estimators. For the robust estimator it also returns the estimate of the residual standard error.
#'
#' @param object an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.
#' @param ... additional other arguments.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
#' @aliases summary.margint summary.margint.cl summary.margint.rob
summary.margint <- function(object,...){
  NextMethod()
}


#' @export
summary.margint.cl <- function(object,...){
  message("Estimate of the intercept: ", round(object$mu,5))
  res <- residuals(object)
  message("Residuals:")
  summary(res)
}


#' @export
summary.margint.rob <- function(object,...){
  message("Estimate of the intercept: ", round(object$mu,5))
  message("Estimate of the residual standard error: ", round(object$sigma,5))
  res <- residuals(object)
  message("Residuals:")
  summary(res)
}


