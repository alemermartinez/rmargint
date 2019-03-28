#' @useDynLib RMI, .registration = TRUE
#' @import stats graphics


#' @title Derivative of Tukey's bi-square loss function
#' @description Derivative of Tukey's bi-square loss function.
#' @param r A real number.
#' @param k Tuning constant.
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @return A real number.
#' @export

psi.tukey <- function(r, k=4.685){
  u <- abs(r/k)
  w <- r*((1-u)*(1+u))^2
  w[u>1] <- 0
  return(w)
}



#' @title Tukey bi-square weight function
#' @description Tukey bi-square weight function.
#' @param r A real number.
#' @param k Tuning constant.
#' @details he weight function used in the re-weighted least squares iterations.
#' @return A real number.
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @export

#Tukey's weight function "Psi(r)/r"
psi.w <- function(r, k= 4.685){
  u <- abs(r/k)
  w <- ((1 + u) * (1 - u))^2
  w[u > 1] <- 0
  return(w)
}


#' @title Derivative of Huber's loss function
#' @description Derivative of Huber's loss function.
#' @param r A real number.
#' @param k Tuning constant.
#' @details The weight function used in the re-weighted least squares iterations.
#' @return A real number.
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @export

# Huber's Psi
psi.huber <- function(r, k=1.345)
  pmin(k, pmax(-k, r))


#' @title Huber weight function
#' @description Huber weight function.
#' @param r A real number.
#' @param k Tuning constant.
#' @details The weight function used in the re-weighted least squares iterations.
#' @return A real number.
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @export


#Huber's weight function "Psi(r)/r"
psi.huber.w <- function(r, k=1.345)
  pmin(1, k/abs(r))


#' @title Euclidean norm of a vector
#' @description Euclidean norm of a vector.
#' @param x A real vector.
#' @return The Euclidean norm of the input vector.
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @export


#Euclidean norm
my.norm.2 <- function(x) sqrt(sum(x^2))


#' @title Epanechnikov kernel
#' @description Epanechnikov kernel.
#' @usage k.epan(x)
#' kernel2(x)
#' @param x A real number.
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @return 0 if abs(x) > 1 and 0.75 * (1 - x^2) otherwise.
#' @export

#Epanechnikov kernel
k.epan<-function(x) {
  a <- 0.75*(1-x^2)
  tmp <- a*(abs(x)<=1)
  return(tmp)
}


#' @export

#Order 2 kernel = Epanechnikov kernel
kernel2<-function(t){
  nucleo<- (3/4)*(1-t^2)*(abs(t)<=1)
  nucleo
}

#- Higher order kernels -#


#' @title Order 4 kernel
#' @description A kernel of order 4.
#' @param x A real number.
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @return 0 if abs(x) > 1 and ( 15/32 ) * ( 1 - x^2 ) * ( 3 - 7 * x^2 ) otherwise.
#' @details A kernel L is a kernel of order 4 if it integrates 1, the integrals of u^j L(u) are 0 for 1 <= j < 4 (j integer) and the integral of u^4 L(u) is different from 0.
#' @export

#Order 4
kernel4<-function(x) {
  a <- (15/32)*(1-x^2)*(3-7*x^2)
  tmp <- a*(abs(x)<= 1)
  return(tmp)
}


#' @title Order 6 kernel
#' @description A kernel of order 6.
#' @param x A real number.
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @return 0 if abs(x) > 1 and ( 105/256 ) * ( 1 - x^2 ) * ( 5 - 30 * x^2 + 33 * x^4 ) otherwise.
#' @details A kernel L is a kernel of order 6 if it integrates 1, the integrals of u^j L(u) are 0 for 1 <= j < 6 (j integer) and the integral of u^6 L(u) is different from 0.
#' @export

#Order 6
kernel6<-function(x) {
  a <- (105/256)*(1-x^2)*(5-30*x^2+33*x^4)
  tmp <- a*(abs(x)<=1)
  return(tmp)
}

#' @title Order 8 kernel
#' @description A kernel of order 8.
#' @param x A real number.
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @return 0 if abs(x) > 1 and ( 315/4096 ) * ( 1 - x^2 ) * ( 35 - 385 * x^2 + 1001 * x^4 - 715 * x^6 ) otherwise.
#' @details A kernel L is a kernel of order 8 if it integrates 1, the integrals of u^j L(u) are 0 for 1 <= j < 8 (j integer) and the integral of u^8 L(u) is different from 0.
#' @export


#Order 8
kernel8<-function(x) {
  a <- (315/4096)*(1-x^2)*(35-385*x^2+1001*x^4-715*x^6)
  tmp <- a*(abs(x)<=1)
  return(tmp)
}

#' @title Order 10 kernel
#' @description A kernel of order 10.
#' @param x A real number.
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @return 0 if abs(x) > 1 and 0.75 * ( 1 - x^2 ) * ( 315/128 - 105/32 * x^2 + 63/64 * x^4 - 3/32 * x^6 - 1/384 * x^8 ) otherwise.
#' @details A kernel L is a kernel of order 10 if it integrates 1, the integrals of u^j L(u) are 0 for 1 <= j < 10 (j integer) and the integral of u^10 L(u) is different from 0.
#' @export

#Order 10
kernel10<-function(x) {
  a <- 0.75*(1-x^2)*(315/128-105/32*x^2+63/64*x^4-3/32*x^6-1/384*x^8)
  tmp <- a*(abs(x)<=1)
  return(tmp)
}

#' @title Classic marginal integration procedures for additive models
#' @description Standard marginal integration procedures for additive models.
#' @param Xp Matrix of explanatory variables (n by p).
#' @param yp  Vector of responses (missing values are allowed).
#' @param point Matrix of points where predictions will be computed and returned.
#' @param windows Vector or a squared matrix of bandwidths for the smoothing estimation procedure.
#' @param epsilon Convergence criterion.
#' @param prob Probabilities of observing each response (n). Defaults to ``NULL''.
#' @param type Three different type of estimators can be selected: type '0' (local constant on all the covariates), type '1' (local linear smoother on all the covariates), type 'alpha' (local polynomial smoother only on the direction of interest).
#' @param degree Degree of the local polynomial smoother in the direction of interest when using the estimator of type 'alpha'. Defaults to ``NULL'' for the case when using estimators of type '0' or '1'.
#' @param orderkernel Order of the kernel used in the nuisance directions when using the estimator of type 'alpha'. Defaults to ``2''.
#' @param qderivate If TRUE, it calculates g^(q+1)/(q+1)! for each component only for the type 'alpha' method. Defaults to ``FALSE''.
#' @param Qmeasure A matrix of points where the integration procedure ocurrs. Defaults to ``NULL'' for calcuting the integrals over the sample.
#' @details Three types of classical marginal integration procedures for additive models, that is, considering a squared loss function.
#' @return
#' \item{mu}{Estimate for the intercept.}
#' \item{g.matrix}{Matrix of estimated additive components (n by p).}
#' \item{prediction }{Matrix of estimated additive components for the points listed in the argument point.}
#' \item{mul}{A vector of size p showing in each component the estimated intercept that considers only that direction of interest when using the type 'alpha' method.}
#' \item{g.derivative }{Matrix of estimated derivatives of the additive components (only when qderivate is ``TRUE'') (n by p).}
#' \item{prediction.derivate }{Matrix of estimated derivatives of the additive components for the points listed in the argument point (only when qderivate is ``TRUE'').}
#' \item{Xp}{Matrix of explanatory variables.}
#' \item{yp}{Vector of responses.}
#'
#' @keywords Additive Models
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @export

## Classic Marginal Integration
margint.cl <- function(Xp, yp, point=NULL, windows, epsilon=1e-6, prob=NULL,
                       type='0', degree=NULL, qderivate=FALSE, orderkernel=2,
                       Qmeasure=NULL) {
  # Xp = covariance matrix (n x q).
  # yp = response vector (NA's are allowed).
  # point = vector of length q or a matrix with q columns where predictions
  #  are computed. If missing, predictions are computed for each row of Xp.
  # windows = vector of length q or a q times q matrix with kernel windows.
  #  The matrix is used for the 'alpha' procedure.
  # epsilon = convergence criterion.
  # prob = probabilities of observing each response (n).
  # type= '0' (local constant), '1' (local linear), 'alpha' (local
  #  polynomial using 'degree').
  # degree = degree of the local polynomial smoother in the direction of
  # interest when using the estimator of type 'alpha'.
  # orderkernel = order of the kernel used in the nuisance directions when
  #  using the estimator of type 'alpha'.
  # qderivate = if TRUE, it calculates g^(q+1)/(q+1)! for each component
  #  only for the type 'alpha' method.
  # Qmeasure = if NULL, the integration procedure is over the sample.

  if(!is.null(dim(Xp))){
    if(type=='alpha'){
      if(is.null(degree)){
        stop("Degree of local polynomial missing")
      }else{
        if( length(windows)==1 ){
          stop("Windows should be a vector o a matrix")
        }
      }
    }else{
      if( (is.matrix(windows)) | (length(windows)==1)  ){
        stop("Windows should be a vector")
      }
    }
  }else{
    if(!is.null(dim(windows))){
      stop("Windows should be a number")
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
    aa <- rep(0,n)

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

#' @title Robust marginal integration procedures for additive models
#' @description Robust marginal integration procedures for additive models.
#' @param Xp Matrix of explanatory variables (n by p).
#' @param yp  Vector of responses (missing values are allowed).
#' @param point Matrix of points where predictions will be computed and returned.
#' @param windows Vector or a squared matrix of bandwidths for the smoothing estimation procedure.
#' @param prob Probabilities of observing each response (n). Defaults to ``NULL''.
#' @param sigma.hat Estimate of the residual standard error. If NULL we use the mad of the residuals obtained with local medians.
#' @param win.sigma Vector of bandwidths for estimating sigma.hat. If NULL it uses the argument windows if it is a vector or its diagonal if it is a matrix.
#' @param epsilon Convergence criterion.
#' @param type Three different type of estimators can be selected: type '0' (local constant on all the covariates), type '1' (local linear smoother on all the covariates), type 'alpha' (local polynomial smoother only on the direction of interest).
#' @param degree Degree of the local polynomial smoother in the direction of interest when using the estimator of type 'alpha'. Defaults to ``NULL'' for the case when using estimators of type '0' or '1'.
#' @param typePhi One of either ``Tukey'' or ``Huber''.
#' @param k.h Tuning constant for a Huber-type loss function. Defaults to ``1.345''.
#' @param k.t Tuning constant for a Tukey-type loss function. Defaults to ``4.685''.
#' @param max.it Maximum number of iterations for the algorithm.
#' @param qderivate If TRUE, it calculates g^(q+1)/(q+1)! for each component only for the type 'alpha' method. Defaults to ``FALSE''.
#' @param orderkernel Order of the kernel used in the nuisance directions when using the estimator of type 'alpha'. Defaults to ``2''.
#' @param Qmeasure A matrix of points where the integration procedure ocurrs. Defaults to ``NULL'' for calcuting the integrals over the sample.
#' @details Three types of robust marginal integration procedures for additive models.
#' @return
#' \item{mu }{Estimate for the intercept.}
#' \item{g.matrix }{Matrix of estimated additive components (n by p).}
#' \item{sigma.hat }{Estimate of the residual standard error.}
#' \item{prediction }{Matrix of estimated additive components for the points listed in the argument point.}
#' \item{mul }{A vector of size p showing in each component the estimated intercept that considers only that direction of interest when using the type 'alpha' method.}
#' \item{g.derivative }{Matrix of estimated derivatives of the additive components (only when qderivate is ``TRUE'') (n by p).}
#' \item{prediction.derivate }{Matrix of estimated derivatives of the additive components for the points listed in the argument point (only when qderivate is ``TRUE'').}
#' \item{Xp}{Matrix of explanatory variables.}
#' \item{yp}{Vector of responses.}
#'
#' @keywords Additive Models
#' @author Alejandra Martinez, Matias Salibian-Barrera
#' @export



## Robust Marginal Integration
margint.rob <- function(Xp, yp, point=NULL, windows, prob=NULL, sigma.hat=NULL,
                        win.sigma=NULL, epsilon=1e-06, type='0', degree=NULL, typePhi='Huber',
                        k.h=1.345, k.t = 4.685, max.it=20, qderivate=FALSE, orderkernel=2,
                        Qmeasure=NULL){

  # Xp = covariance matrix (n x q).
  # yp = response vector (NA's are allowed).
  # point = vector of length q or a matrix with q columns where predictions
  #  are computed. If missing, predictions are computed for each row of Xp.
  # windows = vector of length q x q times q matrix with kernel windows.
  #  The matrix is used for the 'alpha' procedure.
  # epsilon = convergence criterion.
  # prob = probabilities of observing each response (n).
  # type= '0' (local constant), '1' (local linear), 'alpha' (local
  #  polynomial using 'degree').
  # degree = degree of the local polynomial smoother in the direction of
  # interest when using the estimator of type 'alpha'.
  # orderkernel = order of the kernel used in the nuisance directions when
  #  using the estimator of type 'alpha'.
  # qderivate = if TRUE, it calculates g^(q+1)/(q+1)! for each component (only for the
  #  type 'alpha' method.
  # Qmeasure = if NULL, the integration procedure is over the sample.
  # sigma.hat = estimate of the residual standard error. If missing we use the
  #   mad of the residuals obtained with local medians.
  # initial = initial intercept estimate.
  # k.h = tuning constant for the Huber function.
  # k.t = tuning constant for the Tukey function.
  # typePhi = 'Huber' or 'Tukey'
  # max.it = max number of iterations for the robust estimation step.

  if(!is.null(dim(Xp))){
    if(type=='alpha'){
      if(is.null(degree)){
        stop("Degree of local polynomial missing")
      }else{
        if( length(windows)==1 ){
          stop("Windows should be a vector or a matrix")
        }
      }
    }else{
      if( (is.matrix(windows)) | (length(windows)==1) ){
        stop("Windows should be a vector")
      }
    }
  }else{
    if(!is.null(dim(windows))){
      stop("Windows should be a number")
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


#S3 Methods

#' Residuals of a marginal integration fit in an additive model
#'
#' This function returns the residuals of the fitted additive model using one of the three
#' classical or robust marginal integration estimators, as computed with \code{\link{margint.cl}} or
#' \code{\link{margint.rob}}.
#'
#' @param object an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.
#' @param ... additional other arguments.
#'
#' @return A vector of residuals.
#'
#' @author Alejandra Mercedes Martinez \email{alemartinez@unlu.edu.ar}
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
#' @param ... additional other arguments.
#'
#' @return A vector of fitted values.
#'
#' @author Alejandra Mercedes Martinez \email{alemartinez@unlu.edu.ar}
#'
#' @export

predict.margint <- function(object, ...){
  return( rowSums(object$g.matrix) + object$mu )
}

#' Diagnostic plots for objects of class \code{margint}
#'
#' Plot method for class \code{margint}.
#'
#' @param object an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.
#' @param which vector of indices of explanatory variables for which partial residuals plots will be generated. Defaults to all available explanatory variables.
#' @param ask logical value. If \code{TRUE}, the graphical device will prompt before going to the next page/screen of output.
#' @param ... additional other arguments.
#'
#' @author Alejandra Mercedes Martinez \email{alemartinez@unlu.edu.ar}
#'
#' @export
plot.margint <- function(object, which=1:np, ask=FALSE,...){
  Xp <- object$Xp
  np <- dim(Xp)[2]
  opar <- par(ask=ask)
  on.exit(par(opar))
  these <- rep(FALSE, np)
  these[ which ] <- TRUE
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
      plot(Xp[ord,i],object$g.matrix[ord,i],type="l",lwd=3,main="",xlab=x_name,ylab=y_name, ylim=lim_cl)
      points(Xp[,i], res, pch=20,col='gray45')
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
#' @author Alejandra Mercedes Martinez \email{alemartinez@unlu.edu.ar}
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


