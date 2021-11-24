library(PMA)
library(RGCCA)
library(Rcpp)
library(rifle)
library(parallel)
sourceCpp("~/Dropbox/CrossCorrespondance/mutilpleSetLinear/code/msCCA2/src/solvers.cpp")
source("/Users/lg689/Dropbox/CrossCorrespondance/mutilpleSetLinear/code/msCCA2/R/helpers.R")
source("/Users/lg689/Dropbox/CrossCorrespondance/mutilpleSetLinear/code/msCCA2/R/msCCA.R")
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

p = 300; nte = 1000; ncomp = 3; alpha = 0; D = 4; ncomp1 =10 #ncomp-1
n=300; s = 15;D = 4;type = "spiked";redundant = T; seed = 1;#seed = sample(1:10000,1);  

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



D = length(xlist)
ps = sapply(xlist,function(z) dim(z)[2])
pss = c(0,cumsum(ps))
ptotal = pss[D+1]
n = nrow(xlist[[1]])

A = t(xagg)%*%xagg/n
B = array(0, dim = dim(A))
for(d in 1:D){
  ll = (pss[d]+1):pss[d+1]
  B[ll,ll] = t(xlist[[d]])%*%xlist[[d]]/n
}

ptotal = pss[D+1]
A1 = t(xagg.te)%*%xagg.te/nte
B1 = array(0, dim = dim(A1))
for(d in 1:D){
  ll = (pss[d]+1):pss[d+1]
  B1[ll,ll] = t(xlist.te[[d]])%*%xlist.te[[d]]/nte
}

init_method = "soft-thr"; nfolds = 7; foldid = NULL;
rho_maxit = 1000; print_out = 100; l1proximal_maxit = 1e4; rho_tol = 1e-3;
eps =log(ptotal)/n; mode = "full";trace = T


msCCAproximal_l1 = msCCAl1$new(X = xlist, beta_init =NULL, norm_type = 1, init_method = init_method,
                            rho_maxit = rho_maxit, print_out = print_out,
                            l1proximal_maxit = l1proximal_maxit, rho_tol = rho_tol,
                            eps = eps, trace = trace)
msCCAproximal_l1$eps = 1/n
msCCAproximal_l1$eta_ratio = sqrt(1/n)
msCCAproximal_l1$eta_low =sqrt(1/n)*2
msCCAproximal_l1$eta = 0.05
l1norm_max = sum(sapply(msCCAproximal_l1$beta_init,function(z) sum(abs(z))))
l1norm_min = sqrt(2) #sum(abs(U.truth)[,1])
l1norm_max  = min(sqrt(n/5),sqrt(ptotal)) #max(l1norm_max , min(sqrt(n/log(ptotal)),sqrt(ptotal)))
msCCAproximal_l1$print_out = 100
rho_tol = 1e-3;
msCCAproximal_l1$rho_tol = rho_tol
msCCAproximal_l1$rho_maxit = 5000

###########penalized rhos

#####first direction
out = msCCAproximal_l1$direction_update_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R, 
                                                  beta_init =  msCCAproximal_l1$beta_init, 
                                                  l1norm_max =  l1norm_max,   l1norm_min = l1norm_min,
                                                  warm_up = 50, trace = trace)

out1 = msCCAproximal_l1$direction_selection(method = "penalized objective", penalty.C = 2,nfolds = 10, foldid = NULL, seed = 2021, n.core = NULL, trace = F)

step_idx = which.max( out1$evaluation_obj)


out.fix1 = msCCAproximal_l1$direction_update_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R, 
                                                    beta_init =  msCCAproximal_l1$beta_init, 
                                                    l1norm_max =  out$bounds[step_idx],   l1norm_min = out$bounds[step_idx],
                                                    warm_up = 50, trace = trace, record = F)

msCCAproximal_l1$direction_grow(step_idx=step_idx)

beta_selected = msCCAproximal_l1$prev_directions_agg



print(t(beta_selected)%*%A%*%beta_selected/(t(beta_selected)%*%B%*%beta_selected))
print(t(beta_selected)%*%A1%*%beta_selected/(t(beta_selected)%*%B1%*%beta_selected))

beta.fix1 = out.fix1$beta_augs[,ncol(out.fix1$beta_augs)]
print(t(beta.fix1)%*%A%*%beta.fix1/(t(beta.fix1)%*%B%*%beta.fix1))
print(t(beta.fix1)%*%A1%*%beta.fix1/(t(beta.fix1)%*%B1%*%beta.fix1))


#####second direction


out = msCCAproximal_l1$direction_update_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R, 
                                               beta_init =  msCCAproximal_l1$beta_init, 
                                               l1norm_max =  l1norm_max,   l1norm_min = l1norm_min,
                                               warm_up = 50, trace = trace)

out1 = msCCAproximal_l1$direction_selection(method = "penalized objective", penalty.C = 2,nfolds = 10, foldid = NULL, seed = 2021, n.core = NULL, trace = F)


step_idx = which.max( out1$evaluation_obj)

msCCAproximal_l1$direction_grow(step_idx=step_idx)

beta_selected = msCCAproximal_l1$prev_directions_agg

out.fix1 = msCCAproximal_l1$direction_update_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R, 
                                                    beta_init =  msCCAproximal_l1$beta_init, 
                                                    l1norm_max =  out$bounds[step_idx],   l1norm_min = out$bounds[step_idx],
                                                    warm_up = 50, trace = trace, record = F)


beta.fix1 = out.fix1$beta_augs[,ncol(out.fix1$beta_augs)]

print(t(beta_selected)%*%A%*%beta_selected/(t(beta_selected)%*%B%*%beta_selected))
print(t(beta_selected)%*%A1%*%beta_selected/(t(beta_selected)%*%B1%*%beta_selected))


print(t(beta.fix1)%*%A%*%beta.fix1/(t(beta.fix1)%*%B%*%beta.fix1))
print(t(beta.fix1)%*%A1%*%beta.fix1/(t(beta.fix1)%*%B1%*%beta.fix1))

##############################################################cross validation

msCCAproximal_l1 = msCCAl1$new(X = xlist, beta_init =NULL, norm_type = 1, init_method = init_method,
                               rho_maxit = rho_maxit, print_out = print_out,
                               l1proximal_maxit = l1proximal_maxit, rho_tol = rho_tol,
                               eps = eps, trace = trace)
msCCAproximal_l1$eps = 1/n
msCCAproximal_l1$eta_ratio = sqrt(1/n)
msCCAproximal_l1$eta_low =sqrt(1/n)*2
msCCAproximal_l1$eta = sqrt(1/n)
l1norm_max = sum(sapply(msCCAproximal_l1$beta_init,function(z) sum(abs(z))))
l1norm_min = sqrt(2) #sum(abs(U.truth)[,1])
l1norm_max  = min(sqrt(n/5),sqrt(ptotal)) #max(l1norm_max , min(sqrt(n/log(ptotal)),sqrt(ptotal)))
msCCAproximal_l1$print_out = 100
rho_tol = 1e-3;
msCCAproximal_l1$rho_tol = rho_tol
msCCAproximal_l1$rho_maxit = 5000
nfolds = 16
#####first direction
out = msCCAproximal_l1$direction_update_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R, 
                                               beta_init =  msCCAproximal_l1$beta_init, 
                                               l1norm_max =  l1norm_max,   l1norm_min = l1norm_min,
                                               warm_up = 50, trace = trace)

out1 = msCCAproximal_l1$direction_selection(method = "cv", nfolds = nfolds, foldid = NULL, seed = 2022, n.core = NULL)

step_idx = which.max( out1$evaluation_obj)

Zs = array(NA, dim = c(nte, D, ncol(xagg.te%*%out$beta_augs)))
for(d in 1:D){
  ll = (pss[d]+1):pss[d+1]
  Zs[,d,] = xagg.te[,ll]%*%out$beta_augs[ll,]
}
Zsums = apply(Zs,c(1,3), sum)
rho.te.hat = apply(Zsums^2,2,sum)/apply(Zs^2,3,sum)
plot(out1$evaluation_obj,rho.te.hat)
print(step_idx)


out.fix1 = msCCAproximal_l1$direction_update_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R, 
                                                    beta_init =  msCCAproximal_l1$beta_init, 
                                                    l1norm_max =  out$bounds[step_idx],   l1norm_min = out$bounds[step_idx],
                                                    warm_up = 50, trace = trace, record = F)


msCCAproximal_l1$direction_grow(step_idx=step_idx)

beta_selected = msCCAproximal_l1$prev_directions_agg

beta.fix1 = out.fix1$beta_augs[,ncol(out.fix1$beta_augs)]

print(t(beta_selected)%*%A%*%beta_selected/(t(beta_selected)%*%B%*%beta_selected))
print(t(beta_selected)%*%A1%*%beta_selected/(t(beta_selected)%*%B1%*%beta_selected))


print(t(beta.fix1)%*%A%*%beta.fix1/(t(beta.fix1)%*%B%*%beta.fix1))
print(t(beta.fix1)%*%A1%*%beta.fix1/(t(beta.fix1)%*%B1%*%beta.fix1))


#####second direction

out = msCCAproximal_l1$direction_update_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R, 
                                               beta_init =  msCCAproximal_l1$beta_init, 
                                               l1norm_max =  l1norm_max,   l1norm_min = l1norm_min,
                                               warm_up = 50, trace = trace)

out1 = msCCAproximal_l1$direction_selection(method = "cv", nfolds = nfolds, foldid = NULL, seed = 2022, n.core = NULL)


step_idx = which.max( out1$evaluation_obj)

msCCAproximal_l1$direction_grow(step_idx=step_idx)

beta_selected = msCCAproximal_l1$prev_directions_agg

out.fix1 = msCCAproximal_l1$direction_update_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R, 
                                                    beta_init =  msCCAproximal_l1$beta_init, 
                                                    l1norm_max =  out$bounds[step_idx],   l1norm_min = out$bounds[step_idx],
                                                    warm_up = 50, trace = trace, record = F)


beta.fix1 = out.fix1$beta_augs[,ncol(out.fix1$beta_augs)]

print(t(beta_selected)%*%A%*%beta_selected/(t(beta_selected)%*%B%*%beta_selected))
print(t(beta_selected)%*%A1%*%beta_selected/(t(beta_selected)%*%B1%*%beta_selected))


print(t(beta.fix1)%*%A%*%beta.fix1/(t(beta.fix1)%*%B%*%beta.fix1))
print(t(beta.fix1)%*%A1%*%beta.fix1/(t(beta.fix1)%*%B1%*%beta.fix1))


####################################
beta_rifle_init = c()
for(d in 1:D){
  beta_rifle_init  = c(beta_rifle_init, msCCAproximal_l1$beta_init[[d]])
}


tmp1 = rifle(A =A , B = B, init = beta_rifle_init,k = s*D)

print(t(tmp1)%*%A%*%tmp1/(t(tmp1)%*%B%*%tmp1))



print(t(tmp1)%*%A1%*%tmp1/(t(tmp1)%*%B1%*%tmp1))
