require(MASS)
require(stats)
require(geigen)
require(pracma)
Soft <-
  function(a,b){
    if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
    return(sign(a)*pmax(0,abs(a)-b))
  }

updatePi <-
  function(B,sqB,A,H,Gamma,nu,rho,Pi,tau){
    C <- Pi + 1/tau*A-nu/tau*B%*%Pi%*%B+nu/tau*sqB%*%(H-Gamma)%*%sqB
    D <- rho/tau
    return(Soft(C,D))
  }

updateH <-
  function(sqB,Gamma,nu,Pi,K){
    
    temp <- 1/nu * Gamma + sqB%*%Pi%*%sqB
    temp <- (temp+t(temp))/2
    svdtemp <- eigen(temp)
    d <- svdtemp$values
    p <- length(d)
    if(sum(pmin(1,pmax(d,0)))<=K){
      dfinal <- pmin(1,pmax(d,0))
      return(svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors))
    }
    fr <- function(x){
      sum(pmin(1,pmax(d-x,0)))
    }
    # Vincent Vu Fantope Projection
    knots <- unique(c((d-1),d))
    knots <- sort(knots,decreasing=TRUE)
    temp <- which(sapply(knots,fr)<=K)
    lentemp <- tail(temp,1)
    a=knots[lentemp]
    b=knots[lentemp+1]
    fa <- sum(pmin(pmax(d-a,0),1))
    fb <- sum(pmin(pmax(d-b,0),1))
    theta <- a+ (b-a)*(K-fa)/(fb-fa)
    dfinal <- pmin(1,pmax(d-theta,0))
    res <- svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors)
    return(res)
  }

hard <-
  function(U, k){
    if(k>1){
      truncate.value <- sort(apply(U, 1, FUN = function(x) sum(x^2)),decreasing=TRUE)[k]
      U[which(apply(U, 1, FUN = function(x) sum(x^2))<truncate.value)] <- 0
    }else{
      truncate.value <- sort(abs(U),decreasing=TRUE)[k]
      U[which(abs(U)<truncate.value)] <- 0
    }
    return(U)
  }

#'@export
#'\examples{
#'}
subdistance <- function(A, B){
  svdresult = svd(t(A) %*% B);
  U = svdresult$u
  V = svdresult$v
  O = U %*% t(V);
  l = norm(A %*% O-B, type = c('F'));
  return(l)
}

# function to convert Pi from Fantope to input of TGD
# Inputs:
# =======
# Pi:         Output of sgca_init$Pi
# r:          Latent dimension

# Outputs:
# ========
# ainit:   Initialization for the generalized eigenspace
init_process <-
  function(Pi, r){
    ainit = svd(Pi)
    uest <- ainit$u
    dest <- diag(ainit$d)
    if (r == 1){
      ainit <- uest[,1] * sqrt(dest[1:r,1:r])
    } else
      ainit <- uest[,1:r] %*% sqrtm(dest[1:r,1:r])$B
    return(ainit)
  }


#' sgca initialization Inputs:
#' =======
#' A, B:       Pair of matrix to calculate generalized eigenspace (sigma, sigma0)
#' nu:         Parameter of ADMM, default set to 1
#' K:          nuclear norm constraint, equal to r
#' rho:     penalty parameter on the l_1 norm of the solution, scaled by
#             sqrt(log(max(p1,p2))/n)
#' epsilon:    tolerance level for convergence in ADMM
#' maxiter:    maximum number of iterations in ADMM
#' trace:      if set to True will print all iterations 

#' Outputs:
# ========
# $Pi:     optimum of the convex program
#'
#'@export
#'\examples{
#'}
sgca_init <-
  function(A,B,rho,K,nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE){
    p <- nrow(B)
    eigenB <- eigen(B)
    sqB <- eigenB$vectors%*%sqrt(diag(pmax(eigenB$values,0)))%*%t(eigenB$vectors)	
    tau <- 4*nu*eigenB$values[1]^2	
    criteria <- 1e10
    i <- 1
    # Initialize parameters
    H <- Pi <- oldPi <-  diag(1,p,p)
    Gamma <- matrix(0,p,p)
    # While loop for the iterations
    while(criteria > epsilon && i <= maxiter){
      #for (ii in 1:20){
        #Pi <- updatePi(B,sqB,A,H,Gamma,nu,rho,Pi,tau)
      #}
      Pi <- updatePi(B,sqB,A,H,Gamma,nu,rho,Pi,tau)
      
      H <- updateH(sqB,Gamma,nu,Pi,K)
      Gamma <- Gamma + (sqB%*%Pi%*%sqB-H) * nu	
      criteria <- sqrt(sum((Pi-oldPi)^2))
      oldPi <- Pi
      i <- i+1
      if(trace==TRUE)
      {
        print(i)
        print(criteria)
      }
    }
    return(list(Pi=Pi,H=H,Gamma=Gamma,iteration=i,convergence=criteria))
    
  }

#' Function for Thredholded Gradient Descent
#' Inputs:
#' =======
#' A, B:         Pair of matrix to calculate generalized eigenspace (sigma, sigma0)
#' r:            Latent dimension
#' k:            sparsity
#' init:         Initialize estimator of generlized eigenspace, obtained from sgca_int. Need to be k-sparse
#' lambda:       penalty parameter of Lagrangian function f(L), default to be 0.01
#' convergence:  tolerance level for convergence in TGD
#' maxiter:      maximum number of iterations in TGD
#' plot:         if set to True will plot intermediate iterations, need to specify scale variable V (as shown in paper)
#'               default set to False

#' Outputs:
#' ========
#' final_estimate:   final estimation of leading r sparse generalized eigenspace
#'
#'@export
#'\examples{
#'}
sgca_tgd <-
  function(A, B, r, init, k, lambda = 0.01, eta=0.01, convergence=1e-3, maxiter=10000, plot = FALSE){
    #perform hard thresholding
    init <- hard(init, k)
    u <- init
    criteria <- 1e10
    iter <- 0
    error <- rep(0, maxiter)
    # renormalization 
    ut <- init %*% sqrtm(diag(r)+t(u) %*% A %*% u/lambda)$B;
    
    while(criteria > convergence && iter <= maxiter){
      #print(paste0(iter, ":", criteria))
      #perform gradient descent
      grad <- -A %*% ut + lambda * B %*% ut %*% (t(ut) %*% B %*% ut- diag(r));
      vt <- ut - eta * grad
      
      
      # Perform truncation
      vt <- hard(vt, k)
      
      criteria <- sqrt(sum((ut-vt)^2))
      
      ut <- vt
      iter <- iter+1
      if (plot){
        error[iter] <- subdistance(vt, scale)
      }
    }
    if (plot){
      plot(error[1:iter], type='l',  main="Matrix distance and iterations", 
           xlab="Number of iterations", ylab="Matrix distance",)
    }
    final_estimate <- ut %*% sqrtm(t(ut) %*% B %*% ut)$Binv
    return(final_estimate)
  }

#' sgcaTGD_single
#' A: Sigma
#' B: Sigma0
#' r: rank
#' k: row sparsity levels
#'@export
sgcaTGD_single = function(xlist, xagg, r, k,lambda = 0.01, eta=0.01, convergence=1e-3, maxiter=10000, plot = FALSE){
  D = length(xlist)
  ps = sapply(xlist, function(z) dim(z)[2])
  pss = c(0,cumsum(ps))
  ptotal = pss[D+1]
  Sigmahat_tmp =t(xagg)%*%xagg/n 
  Lambdahat_tmp = matrix(0, nrow = ptotal, ncol = ptotal)
  for(d in 1:D){
    ll = (pss[d]+1):(pss[d+1])
    Lambdahat_tmp[ll,ll] = t(xlist[[d]])%*%xlist[[d]]/n
  }
  ag = initial.convex(A=Sigmahat_tmp, B=Lambdahat_tmp, lambda = 0.5*sqrt(log(p)/n), K = r, nu = 1, epsilon = 0.05, maxiter = 1000, trace = F)
  #ag = sgca_init(A=Sigmahat_tmp, B=Lambdahat_tmp, rho=0.5*sqrt(log(p)/n),K=r ,nu=1,trace=FALSE)
  ainit <- init_process(ag$Pi, r)
  #ainit = NULL; for(l in 1:4){ainit = rbind(ainit,fitted$fitted_model$prev_directions[[l]])}
  if(is.null(dim(ainit))){
    ainit = matrix(ainit, ncol = 1)
  }
  print('Finish convex initialization.')
  final <- sgca_tgd(A=Sigmahat_tmp, B=Lambdahat_tmp, r=r, init = ainit,k = k,lambda =lambda, eta=eta,convergence=convergence,maxiter= maxiter, plot =FALSE)
  print('Finish final estimation.')
  return(list(ainit = ainit, final = final))
}
#' sgcaTGD_wrapper
#' A: Sigma
#' B: Sigma0
#' r: rank
#' ks: row sparsity levels
#'@export
sgcaTGD_wrapper <- function(xlist, xagg, r, foldid,ks = NULL,
                            lambda = 0.01, eta=0.01, convergence=1e-3, maxiter=10000, plot = FALSE){
  if(is.null(ks)){
    ks = c(0:14) * 5
    ks[1] = 1
  }
  nfolds = max(foldid)
  n = nrow(xlist[[1]])
  D = length(xlist)
  ps = sapply(xlist, function(z) dim(z)[2])
  pss = c(0,cumsum(ps))
  ptotal = pss[D+1]
  # Cross-validation
  cv_denominator = list();
  cv_numerator = list()
  for(ik in 1:length(ks)){
    cv_denominator[[ik]] = matrix(0, r, r)
    cv_numerator[[ik]] = matrix(0, r, r)
  }
  eta_use = eta
  for(fold_id in 1:nfolds){
    print(paste0("start fold ",fold_id))
    train_id = which(foldid != fold_id)
    xlist_tmp = list()
    xagg_temp = xagg[train_id,]
    for(i in 1:D){
      xlist_tmp[[i]] = xlist[[i]][train_id,]
    }
    Sigmahat_tmp =t(xagg_temp)%*%xagg_temp/n 
    Lambdahat_tmp = matrix(0, nrow = ptotal, ncol = ptotal)
    for(d in 1:D){
      ll = (pss[d]+1):(pss[d+1])
      Lambdahat_tmp[ll,ll] = t(xlist_tmp[[d]])%*%xlist_tmp[[d]]/n
    }
    xlist_tmp1 = list()
    test_id = which(foldid == fold_id)
    xagg_temp1 = xagg[test_id,]
    for(i in 1:D){
      xlist_tmp1[[i]] = xlist[[i]][test_id,]
    }
    Sigmahat_tmp1 =t(xagg_temp1)%*%xagg_temp1/n 
    Lambdahat_tmp1 = matrix(0, nrow = ptotal, ncol = ptotal)
    for(d in 1:D){
      ll = (pss[d]+1):(pss[d+1])
      Lambdahat_tmp1[ll,ll] = t(xlist_tmp1[[d]])%*%xlist_tmp1[[d]]/n
    }
    ##########
    syst1 = Sys.time()
    ag = initial.convex(A=Sigmahat_tmp, B=Lambdahat_tmp, lambda = 0.5*sqrt(log(p)/n), K = r, nu = 1, epsilon = 0.05, maxiter = 1000, trace = T)
    #ag = sgca_init(A=Sigmahat_tmp, B=Lambdahat_tmp, rho=0.5*sqrt(log(p)/n),K=r ,nu=1,trace=T)
    ainit <- init_process(ag$Pi, r)
    if(is.null(dim(ainit))){
      ainit = matrix(ainit, ncol = 1)
    }
    syst2 = Sys.time()
    print(paste0("Fold",fold_id,' Finish convex initialization in ', (syst2-syst1), " min"))
    for(ik in 1:length(ks)){
      k = ks[ik]
      syst1 = Sys.time()
      out <- try(sgca_tgd(A=Sigmahat_tmp, B=Lambdahat_tmp,r=r,init=ainit,k=k,lambda =lambda, eta=eta_use,convergence=convergence,maxiter= maxiter, plot =  FALSE))
      while("try-error" %in% class(out)){
        eta_use = eta_use/2
        out <- try(sgca_tgd(A=Sigmahat_tmp, B=Lambdahat_tmp,r=r,init=ainit,k=k,lambda =lambda, eta=eta_use,convergence=convergence,maxiter= maxiter, plot =  FALSE))
      }
      syst2 = Sys.time()
      print(paste0('Finish Fold ', fold_id, ", sparsity=",k, ", in ", (syst2-syst1), " min"))
      denominator = t(out)%*%Lambdahat_tmp1%*%out
      numerator = t(out)%*%Sigmahat_tmp1%*%out
      cv_denominator[[ik]] =cv_denominator[[ik]]+denominator
      cv_numerator[[ik]] =cv_numerator[[ik]]+numerator
    }
    print(paste0("end fold ",fold_id))
  }
  #cv_denominator = cv_denominator+1E-4
  cv_criteria = sapply(cv_numerator, function(z) sum(diag(z)))
  ik0 = which.max(cv_criteria)
  k0 = ks[ik0]
  syst1 = Sys.time()
  final_out <- try(sgcaTGD_single(xlist=xlist, xagg = xagg, r = r, k = k0, 
                          lambda =lambda, eta=eta_use,convergence=convergence,maxiter= maxiter, plot =  plot))
  while("try-error" %in% class(final_out)){
    eta_use = eta_use/2
    final_out <- try(sgcaTGD_single(xlist=xlist, xagg = xagg, r = r, k = k0, 
                                    lambda =lambda, eta=eta_use,convergence=convergence,maxiter= maxiter, plot =  plot))
    
  }
  syst2 = Sys.time()
  print(paste0('Finish final estimation, selected sparsity at',k0, " in ", (syst2-syst1), " min"))
  return(list(final_out=final_out, k = k0, eta_use = eta_use))
}

