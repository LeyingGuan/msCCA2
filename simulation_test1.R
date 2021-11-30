library(PMA)
library(RGCCA)
library(Rcpp)
library(rifle)
library(parallel)
library(msCCA2)
library(parallel)
library(doParallel)
library(doMC)
sourceCpp("~/Dropbox/CrossCorrespondance/mutilpleSetLinear/code/msCCA2/src/solvers.cpp")
source("~/Dropbox/CrossCorrespondance/mutilpleSetLinear/code/msCCA2/R/helpers.R")
source("~/Dropbox/CrossCorrespondance/mutilpleSetLinear/code/msCCA2/R/msCCA.R")


#library(msCCA)
#source("/home/lg689/project/msCCA/algorithm/msCCAfull.R")
#source("/home/lg689/project/msCCA/algorithm/helpers.R")
#source("/Users/lg689/Dropbox/CrossCorrespondance/mutilpleSetLinear/code/msCCAproximal/R/msCCAfull.R")
#source("/Users/lg689/Dropbox/CrossCorrespondance/mutilpleSetLinear/code/msCCAproximal/R/helpers.R")
orthogonalize_func = function(U, Lambda_svd){
  Uorth =(Lambda_svd$u)%*%diag(sqrt(Lambda_svd$d))%*%t(Lambda_svd$u)%*%U
  if(ncol(U)>1){
    for(k in 1:(ncol(U)-1)){
      for(k1 in (k+1):ncol(U)){
        Uorth[,k1] = lm(Uorth[,k1]~Uorth[,k]-1)$residuals
      }
    }
  }
  
  for(k in 1:ncol(U)){
    Uorth[,k] = Uorth[,k]/sqrt(sum(Uorth[,k]^2))
  }
  Utransformed =Lambda_svd$u%*%diag(1.0/sqrt(Lambda_svd$d))%*%t((Lambda_svd$u))%*% Uorth
  return(Utransformed)
}
print(detectCores())
nfolds = min(detectCores()-1,10)
sim_data_mCCA = function(n = 200, nte = 1000, p = 200, s = 10, D = 5, seed =2021, ncomp = 5,
                         redundant =T, type = "identity"){
  set.seed(seed)
  Lambda_list = list()
  ps = rep(p, D)
  pss = c(0, cumsum(ps))
  ptotal = p*D
  for(d in 1:D){
    if(type == "identity"){
      Lambda_list[[d]] = diag(rep(1, p))
    }else if(type == "toplitz"){
      A1 = matrix(rep(1:p, p), byrow = F, ncol = p)
      A2 = matrix(rep(1:p, p), byrow = T, ncol = p)
      A3 = abs(A1-A2)
      Lambda_list[[d]] = 0.3^(A3)
    }else if(type == "sparseInv"){
      
    }else if(type == "spiked"){
      U = matrix(rnorm(p*ncomp),ncol = ncomp)
      U = apply(U,2,function(z) z/sqrt(sum(z^2)))
      Lambda_list[[d]] = diag(rep(1, p))+100*U%*%t(U)
    }else{
      stop("unsupported covariance type.")
    }
  }
  Lambda = matrix(0, ncol = ptotal, nrow = ptotal)
  for(d in 1:D){
    idx = (pss[d]+1):pss[d+1]
    Lambda[idx, idx] = Lambda_list[[d]]
  }
  Lambda_svd = svd(Lambda)
  ##generate U, V: signal variance s * log(p)/n
  U = matrix(rnorm(ptotal*ncomp),ncol = ncomp) 
  ##randomly select some non-zeros from each assays
  for(d in 1:D){
    temp_U =  U[(pss[d]+1):(pss[d+1]),,drop = F]
    for(k in 1:ncomp){
      idx = sample(1:p, s)
      temp_U[-idx,k] = 0 
    }
    temp_U = scale(temp_U, center = F)
    temp_U = temp_U 
    U[(pss[d]+1):(pss[d+1]),] = temp_U
  }
  ###normalizing with respect to Lambda
  for(d in 1:D){
    temp_U =  U[(pss[d]+1):(pss[d+1]),,drop = F]
    if((d<=2)|(d>2 & !(redundant))){
      Lambda_svd1 = svd(Lambda_list[[d]])
      temp_U = orthogonalize_func(temp_U,Lambda_svd1)
      U[(pss[d]+1):(pss[d+1]),] = temp_U
    }else{
      U[(pss[d]+1):(pss[d+1]),] = 0
    }
  }
  if(!redundant){
    rhos = (.9 - c(0:(ncomp-1))/5)
  }else{
    rhos = (.9 - c(0:(ncomp-1))/5)
  }
  if(ncol(U)>1){
    Sigma = (Lambda%*%U)%*%diag(rhos)%*%(t(U)%*%Lambda)
  }else{
    Sigma = matrix((Lambda%*%U),ncol=1)%*%(t(U)%*%Lambda)*rhos
  }
  
  for(d in 1:D){
    idx = (pss[d]+1):pss[d+1]
    Sigma[idx, idx] = Lambda_list[[d]]
  }
  tmp = eigen(Sigma, symmetric = TRUE, only.values = FALSE, EISPACK = FALSE)
  if(min(tmp$values)<0){
    stop("none PD covariance.")
  }
  tmp1 = Lambda_svd$u%*%diag(sqrt(1.0/Lambda_svd$d))%*%t(Lambda_svd$u)
  tmp2 = tmp1%*%Sigma%*%t(tmp1)
  tmp3 =eigen(tmp2, symmetric = T)
  U=tmp3$vectors[,1:ncomp]
  rhos =tmp3$values[1:ncomp]
  U = tmp1%*%U
  Xagg = matrix(rnorm(n=n*ptotal), ncol = ptotal)%*%tmp$vectors%*%diag(sqrt(tmp$values))%*%t(tmp$vectors)
  Xtestagg = matrix(rnorm(n=nte*ptotal), ncol = ptotal)%*%tmp$vectors%*%diag(sqrt(tmp$values))%*%t(tmp$vectors)
  Xagg = scale(Xagg)
  Xtestagg = scale(Xtestagg)
  ##check
  X = list()
  Xte = list()
  Zsum = array(0, dim = c(n,ncomp))
  Zsum.te =  array(0, dim = c(nte,ncomp))
  Z = array(0, dim = c(n,D,ncomp))
  Z.te = array(0, dim = c(nte, D, ncomp))
  for(d in 1:D){
    ll = (pss[d]+1):pss[d+1]
    X[[d]] = Xagg[,ll]
    Xte[[d]] = Xtestagg[,ll]
    Z[,d,] = X[[d]]%*%U[ll,,drop = F]
    Z.te[,d,] = Xte[[d]]%*%U[ll,]
  }
  Zsum = apply(Z,c(1,3),sum)
  Zsum.te = apply(Z.te,c(1,3),sum)
  return(list(X = X, X.agg = Xagg, Xte = Xte, Xte.agg= Xtestagg,
              U = U, rhos = rhos, Z = Z, Z.te = Z.te, Zsum = Zsum, Zsum.te = Zsum.te,
              Sigma = Sigma))
}
p = 300; nte = 1000; ncomp = 3; alpha = 0; D = 4; ncomp1 =ncomp-1
n=500; s = 5;D = 4;type = "identity";redundant = T; seed = 219;#seed = sample(1:10000,1);  

dat = sim_data_mCCA(n = n, nte = nte, p =p, s = s, D = D, seed =seed, ncomp = ncomp,
                    redundant = redundant, type = type)

xlist = dat$X
xlist.te = dat$Xte
xagg = dat$X.agg
xagg.te = dat$Xte.agg
Zsum.truth = dat$Zsum
Zsum.te.truth = dat$Zsum.te
U.truth = dat$U
rhos.truth = dat$rhos
rhos.truth
Sigma =  dat$Sigma
start_time =c()
end_time = c()
nte = dim(xagg.te)[1]

# D = length(xlist)
# ps = sapply(xlist,function(z) dim(z)[2])
# pss = c(0,cumsum(ps))
# ptotal = pss[D+1]
# n = nrow(xlist[[1]])
# 
# A = t(xagg)%*%xagg/n
# B = array(0, dim = dim(A))
# for(d in 1:D){
#   ll = (pss[d]+1):pss[d+1]
#   B[ll,ll] = t(xlist[[d]])%*%xlist[[d]]/n
# }
# 
# 
# A1 = t(xagg.te)%*%xagg.te/nte
# B1 = array(0, dim = dim(A1))
# for(d in 1:D){
#   ll = (pss[d]+1):pss[d+1]
#   B1[ll,ll] = t(xlist.te[[d]])%*%xlist.te[[d]]/nte
# }

#n.core = strtoi(Sys.getenv("SLURM_CPUS_PER_TASK",unset=1))
n.core = detectCores()-1
print(n.core)



multi.core= "doparallel"

nfolds = 5
set.seed(seed)
foldid = sample(rep(1:nfolds, each = ceiling(n/nfolds)), n)
eta = 0.05; maxit = 5000; s_upper = n/4; eta_ratio = 0.05; penalty.C = 1.9;
start_times =rep(NA, 3)
end_times =rep(NA, 3)
start_times[1] = Sys.time()

fitted1 = msCCAl1func(xlist = xlist, ncomp=ncomp1, xlist.te =xlist.te, init_method = "soft-thr", foldid = foldid, penalty.C=2,
                      l1norm_max =sqrt(s_upper), l1norm_min = sqrt(2), eta = eta, eta_ratio = eta_ratio,
                      rho_maxit = maxit, print_out = 100, step_selection = "penalized", seed = 2021, multi.core=multi.core)
end_times[1] = Sys.time()  
print(end_times[1]-start_times[1])
print(fitted1$errors_track_selected)
 
if(multi.core=="doparallel"){
  registerDoParallel(n.core)
}
start_times[2] = Sys.time()

fitted2 = msCCAl1func(xlist = xlist, ncomp=ncomp1, xlist.te =xlist.te, init_method = "soft-thr", foldid = foldid, penalty.C=2,
                      l1norm_max =sqrt( s_upper), l1norm_min = sqrt(2), eta = eta, eta_ratio = eta_ratio,
                      rho_maxit = maxit, print_out = 100, step_selection = "cv", seed = 2021, multi.core=multi.core)
end_times[2] = Sys.time()
print(fitted2$errors_track_selected)
print(end_times[2]-start_times[2])

start_times[3] = Sys.time()
fitted3 = riffle_sequential(xlist = xlist, ncomp = ncomp1, xlist.te = xlist.te, foldid = foldid, maxiter =maxit, eta = eta,
                            ss = floor(seq(sqrt(2), sqrt(s_upper), length.out = 10)^2),  n.core = NULL, seed = seed, multi.core=multi.core)
end_times[3] = Sys.time()  
print(fitted3$errors_track)
print(end_times[3]-start_times[3])

print(end_times-start_times)
print(fitted1$errors_track_selected)
print(fitted2$errors_track_selected)
print(fitted3$errors_track)

cor(fitted1$fitted_model$prev_directions_agg,U.truth)
cor(fitted2$fitted_model$prev_directions_agg,U.truth)
cor(fitted3$betas,U.truth)



