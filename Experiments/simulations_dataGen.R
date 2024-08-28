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



p = 500; nte = 2000; ncomp = 3; alpha = 0; D = 4; ncomp1 = ncomp-1
path_name = "data/"
set.seed(2024)
for(n in c(300, 1000)){
for(s in c(1, 5, 15)){
for(type in c("identity", "toplitz", "spiked")){
for(redundant in c(T, F)){
print(paste0("simulation_data_n", n, "_s",s,"_",type,"_redundant",redundant))
for(iseed in c(1:100)){
seed = sample(1:100000,1)
save_file_name=paste0(path_name,"simulation_data_n", n, "_s",s,"_",type,"_redundant",redundant,"_iseed",iseed, ".rds")
dat = sim_data_mCCA(n = n, nte = nte, p =p, s = s, D = D, seed =seed, ncomp = ncomp,redundant = redundant, type = type)
saveRDS(dat , file = save_file_name)
}
}
}
}
}
