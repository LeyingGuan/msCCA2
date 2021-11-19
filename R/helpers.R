#' rho_estimation: estimates mCCA coefficients without deflation
#'@export
#'\examples{
#'}
rho_estimation = function(xlist, betas){
  a = 0
  b = 0
  D = length(betas)
  n = nrow(xlist[[1]])
  K = ncol(betas[[1]])
  if(is.null(K)){
    K=1
  }
  Z = array(0, dim = c(n, D, K))
  Zsum = array(0, dim = c(n, K))
  for(d in 1:length(betas)){
      Z[,d,] = xlist[[d]]%*%betas[[d]]/sqrt(n)
    }
  Zsum = apply(Z,c(1,3),sum)
  a = apply(Z^2,3,sum)
  b = apply(Zsum^2, 2,sum)
  rho = b/a
  return(list(Z = Z, Zsum = Zsum, rho = rho, a = a, b = b))
}

#' rho_estimation_deflated: estimates mCCA coefficients withdeflation
#'@export
#'\examples{
#'}
rho_estimation_deflated = function(xlist, betas){
  a = 0
  b = 0
  D = length(betas)
  n = nrow(xlist[[1]])
  K = ncol(betas[[1]])
  if(is.null(K)){
    K=1
  }
  Z = array(0, dim = c(n, D, K))
  Zsum = array(0, dim = c(n, K))
  ps = sapply(xlist, function(z) dim(z)[2])
  pss = c(0, cumsum(ps))
  #concatenated vector
  Zs.agg = array(0, dim = c(n*D, K))
  for(d in 1:length(betas)){
      Z[,d,] = xlist[[d]]%*%betas[[d]]/sqrt(n)
      ll = (n*(d-1)+1):(n*d)
      Zs.agg[ll,] = Z[,d,] 
  }
  #deflation
  if(K>1){
    for(k in 1:(K-1)){
      for(k1 in (k+1):K){
        Zs.agg[,k1] = lm(Zs.agg[,k1]~Zs.agg[,k])$residuals
      }
    }
  }
  for(d in 1:length(betas)){
    ll = (n*(d-1)+1):(n*d)
    Z[,d,] =Zs.agg[ll,]
  }
  Zsum = apply(Z,c(1,3),sum)
  a = apply(Z^2,3,sum)
  b = apply(Zsum^2, 2,sum)
  rho = b/a
  return(list(Z = Z, Zsum = Zsum, rho = rho, a = a, b = b))
}


#' my_sggca_cv: sgcca tuning.
#'@export
my_sggca_cv = function(xlist, sparsities_grids, connection = NULL, nfolds = NULL, foldid = NULL,
                       n.core = NULL){
  n = nrow(xlist[[1]])
  D = length(xlist)
  L = nrow(sparsities_grids)
  ps = sapply(xlist, function(z) dim(z)[2])
  if(is.null(n.core)){
    n.core = min(nfolds, detectCores())
  }
  if(!is.null(foldid)){
    nfolds = length(unique(foldid))
  }else{
    tmp = rep(1:nfolds, each = floor(n/nfolds))
    if(length(tmp)<n){
      tmp = c(c(1:(n-length(tmp))),tmp)
    }
    foldid = sample(1:n, n, replace = F)
    foldid = tmp[foldid]
  }
  if(is.null(connection)){
    connection = 1 - diag(length( xlist))
  }
  func0 = function(input0){
    nf = input0[1]
    silcence = input0[2]
    print(paste0("fold", nf))
    Xtmp = list()
    for(d in 1:D){
      Xtmp[[d]] = xlist[[d]][foldid!=nf,]
    }
    tmp = array(0, dim = c(sum(foldid==nf), D,L))
    for(l in 1:L){
      if(silcence==1){
        print(paste0("param:", l))
      }
      sparsity1 = sparsities_grids[l,]
      result.sgcca.tmp = sgcca(A = Xtmp, C = connection,
                               c1= sparsity1, ncomp = rep(1,length(xlist)),
                               scheme = "horst", verbose = F)
      for(d in 1:D){
        tmp[,d, l] = xlist[[d]][foldid==nf,]%*%result.sgcca.tmp$astar[[d]]
      }
    }
    out = list(Zcollections_nf = tmp)
    return(out)
  }
  inputs = list()
  for(nf in 1:nfolds){
    inputs[[nf]] = c(nf, 0)
    if(nf%%n.core == 1){
      inputs[[nf]][2]=1
    }
  }
  Zcollections = array(0, dim = c(n, D, L))
  out_collections = mclapply( inputs, func0, mc.cores = n.core)
  if(class(out_collections[[1]]) == "try-error"){
    for(i in 1:nfolds){
      inputs[[i]][3] = 1
    }
    out_collections = lapply( inputs, func0)
  }
  for(nf in 1:nfolds){
    Zcollections[foldid==nf, ,] = out_collections[[nf]][[1]]
  }
  #calculate rho
  Zsum = apply(Zcollections,c(1,3),sum)
  a = apply(Zsum^2,2,sum)
  b = apply(Zcollections^2,c(3), sum)
  rhos = a/b
  idx = which.max(rhos)
  sparsity0 = sparsities_grids[idx]
  return(list(beta_penalties = sparsity0, rhos = rhos))
}


#' PMA_wrapper: PMA wrapper
#'@export
PMA_wrapper = function(xlist, xlist.te, ncomp, nperms = 10){
  perm.out <- MultiCCA.permute(xlist, type=rep("standard", length(xlist)), nperms = nperms)
  fitted <- MultiCCA(xlist, type=rep("standard", length(xlist)),penalty=perm.out$bestpenalties, ncomponents=ncomp)
  fitted$prev_directions =  fitted$ws
  #naive
  aa = rho_estimation(xlist,fitted$prev_directions); 
  naive_rho.tr = aa$rho; 
  Zs = aa$Z; Zsum = aa$Zsum;
  naive_rho.te = NULL;Zs.te = NULL; Zsum.te = NULL
  if(!is.null(xlist.te)){
    aa = rho_estimation(xlist.te,fitted$prev_directions); 
    naive_rho.te = aa$rho; 
    Zs.te = aa$Z; Zsum.te = aa$Zsum;
  }
  aa = rho_estimation_deflated(xlist,fitted$prev_directions); 
  deflated_rho.tr = aa$rho; 
  Zsum.deflated =aa$Zsum;
  Zs.deflated = aa$Z
  deflated_rho.te = NULL
  Zsum.te.deflated = NULL
  Zs.te.deflated = NULL
  if(!is.null(xlist.te)){
    aa = rho_estimation_deflated(xlist.te,fitted$prev_directions); 
    deflated_rho.te = aa$rho; 
    Zs.te.deflated = aa$Z;  Zsum.te.deflated = aa$Zsum;
  }
  naive_Zs = list(Zs = Zs, Zsum = Zsum, Zs.te = Zs.te, Zsum.te = Zsum.te,
                  rho.tr = naive_rho.tr, rho.te = naive_rho.te)
  deflated_Zs = list(Zs = Zs.deflated, Zsum = Zsum.deflated, Zs.te = Zs.te.deflated, 
                     Zsum.te = Zsum.te.deflated,
                     rho.tr = deflated_rho.tr, rho.te = deflated_rho.te)
  #calculate Zsum
  return(list(fitted_model = fitted, naive_Zs = naive_Zs,  deflated_Zs =  deflated_Zs))
}

#' rgcca_wrapper: rgcca wrapper
#'@export
rgcca_wrapper = function(xlist,xlist.te, ncomp){
  fitted <- try(rgcca(A = xlist,  tau = "optimal",
                            scheme = "horst", verbose =F, ncomp = rep(ncomp,length(xlist))))
  if("try-error"%in%class(fitted)){
    fitted <- try(rgcca(A = xlist,  tau = "optimal",
                              scheme = "horst", verbose =F, 
                              ncomp = rep(ncomp,length(xlist)),
                              init = "random"))
    
  }
  if(!"try-error"%in%class( fitted)){
    fitted$prev_directions =  fitted$astar
  }
  #naive
  if(!"try-error"%in%class(fitted)){
    aa = rho_estimation(xlist,fitted$prev_directions); 
    naive_rho.tr = aa$rho; 
    Zs = aa$Z; Zsum = aa$Zsum;
    naive_rho.te = NULL;Zs.te = NULL; Zsum.te = NULL
    if(!is.null(xlist.te)){
      aa = rho_estimation(xlist.te,fitted$prev_directions); 
      naive_rho.te = aa$rho; 
      Zs.te = aa$Z; Zsum.te = aa$Zsum;
    }
    aa = rho_estimation_deflated(xlist,fitted$prev_directions); 
    deflated_rho.tr = aa$rho; 
    Zsum.deflated =aa$Zsum;
    Zs.deflated = aa$Z
    deflated_rho.te = NULL
    Zsum.te.deflated = NULL
    Zs.te.deflated = NULL
    if(!is.null(xlist.te)){
      aa = rho_estimation_deflated(xlist.te,fitted$prev_directions); 
      deflated_rho.te = aa$rho; 
      Zs.te.deflated = aa$Z;  Zsum.te.deflated = aa$Zsum;
    }
    naive_Zs = list(Zs = Zs, Zsum = Zsum, Zs.te = Zs.te, Zsum.te = Zsum.te,
                    rho.tr = naive_rho.tr, rho.te = naive_rho.te)
    deflated_Zs = list(Zs = Zs.deflated, Zsum = Zsum.deflated, Zs.te = Zs.te.deflated, 
                       Zsum.te = Zsum.te.deflated,
                       rho.tr = deflated_rho.tr, rho.te = deflated_rho.te)
  }else{
    naive_Zs = NULL
    deflated_Zs = NULL
  }

  #calculate Zsum
  return(list(fitted_model = fitted, naive_Zs = naive_Zs,  deflated_Zs =  deflated_Zs))
  
}

#' sgcca_wrapper: sgcca wrapper
#'@export
sgcca_wrapper = function(xlist,xlist.te, ncomp,mfolds = NULL, foldid = NULL){
  D = length(xlist)
  ps = sapply(xlist, function(z) dim(z)[2])
  pss = c(0,cumsum(ps))
  sparsities_grids = matrix(0, ncol = length(xlist), nrow = 10)
  for(d in 1:D){
    sparsities_grids[,d] = sort(exp(seq(from = log(1.0), to =log(sqrt(ps[d]/2)), length.out = 10)))/sqrt(ps[d])
    aa = sparsities_grids[,d];
    aa[aa>1.0] = 1.0
    aa[aa<(1.0/sqrt(ps[d]))] =(1.0/sqrt(ps[d]))[aa<(1.0/sqrt(ps[d]))]
    sparsities_grids[,d] = aa
  }
  tmp  = my_sggca_cv(xlist = xlist, sparsities_grids = sparsities_grids, connection = NULL, nfolds =nfolds, foldid = foldid, n.core = NULL)
  best_penalties =  tmp$beta_penalties
  
  fitted = sgcca(A = xlist, 
                       c1= best_penalties , ncomp = rep(ncomp,length(xlist)),
                       scheme = "horst", verbose = F)
  if(!"try-error"%in%class( fitted)){
    fitted$prev_directions =  fitted$astar
  }
  #naive
  if(!"try-error"%in%class(fitted)){
    aa = rho_estimation(xlist,fitted$prev_directions); 
    naive_rho.tr = aa$rho; 
    Zs = aa$Z; Zsum = aa$Zsum;
    naive_rho.te = NULL;Zs.te = NULL; Zsum.te = NULL
    if(!is.null(xlist.te)){
      aa = rho_estimation(xlist.te,fitted$prev_directions); 
      naive_rho.te = aa$rho; 
      Zs.te = aa$Z; Zsum.te = aa$Zsum;
    }
    aa = rho_estimation_deflated(xlist,fitted$prev_directions); 
    deflated_rho.tr = aa$rho; 
    Zsum.deflated =aa$Zsum;
    Zs.deflated = aa$Z
    deflated_rho.te = NULL
    Zsum.te.deflated = NULL
    Zs.te.deflated = NULL
    if(!is.null(xlist.te)){
      aa = rho_estimation_deflated(xlist.te,fitted$prev_directions); 
      deflated_rho.te = aa$rho; 
      Zs.te.deflated = aa$Z;  Zsum.te.deflated = aa$Zsum;
    }
    naive_Zs = list(Zs = Zs, Zsum = Zsum, Zs.te = Zs.te, Zsum.te = Zsum.te,
                    rho.tr = naive_rho.tr, rho.te = naive_rho.te)
    deflated_Zs = list(Zs = Zs.deflated, Zsum = Zsum.deflated, Zs.te = Zs.te.deflated, 
                       Zsum.te = Zsum.te.deflated,
                       rho.tr = deflated_rho.tr, rho.te = deflated_rho.te)
  }else{
    naive_Zs = NULL
    deflated_Zs = NULL
  }
  
  #calculate Zsum
  return(list(fitted_model = fitted, naive_Zs = naive_Zs,  deflated_Zs =  deflated_Zs))
  
  
}


#' riffle_wrapper: riffle wrapper
#'@export
riffle_wrapper = function(xlist, xagg, xlist.te, msCCAl0_fitted){
  D = length(xlist)
  ps = sapply(xlist, function(z) dim(z)[2])
  pss = c(0,cumsum(ps))
  beta_rifle_init = c()
  riffle_norm = 0
  for(d in 1:D){
    beta_rifle_init  = c(beta_rifle_init, msCCAl0_fitted$beta_inits[[1]][[d]])
  }
  riffle_norm = msCCAl0_fitted$selected_penalties[1]
  A = t(xagg)%*%xagg/n
  B = array(0, dim = dim(A))
  for(d in 1:D){
    ll = (pss[d]+1):pss[d+1]
    B[ll,ll] = t(xlist[[d]])%*%xlist[[d]]/n
  }
  ptotal = pss[D+1]
  tmp1 = rifle(A =A , B = B, init = beta_rifle_init,k = riffle_norm)
  rifle_beta_list = list()
  for(d in 1:D){
    ll = (pss[d]+1):pss[d+1]
    rifle_beta_list[[d]]= tmp1[ll]
  }
  #naive
  aa = rho_estimation(xlist,rifle_beta_list); 
  naive_rho.tr = aa$rho; 
  Zs = aa$Z; Zsum = aa$Zsum;
  naive_rho.te = NULL;Zs.te = NULL; Zsum.te = NULL
  if(!is.null(xlist.te)){
    aa = rho_estimation(xlist.te,rifle_beta_list); 
    naive_rho.te = aa$rho; 
    Zs.te = aa$Z; Zsum.te = aa$Zsum;
  }
  aa = rho_estimation_deflated(xlist,rifle_beta_list); 
  deflated_rho.tr = aa$rho; 
  Zsum.deflated =aa$Zsum;
  Zs.deflated = aa$Z
  deflated_rho.te = NULL
  Zsum.te.deflated = NULL
  Zs.te.deflated = NULL
  if(!is.null(xlist.te)){
    aa = rho_estimation_deflated(xlist.te,rifle_beta_list); 
    deflated_rho.te = aa$rho; 
    Zs.te.deflated = aa$Z;  Zsum.te.deflated = aa$Zsum;
  }
  naive_Zs = list(Zs = Zs, Zsum = Zsum, Zs.te = Zs.te, Zsum.te = Zsum.te,
                  rho.tr = naive_rho.tr, rho.te = naive_rho.te)
  deflated_Zs = list(Zs = Zs.deflated, Zsum = Zsum.deflated, Zs.te = Zs.te.deflated, 
                     Zsum.te = Zsum.te.deflated,
                     rho.tr = deflated_rho.tr, rho.te = deflated_rho.te)
  #calculate Zsum
  return(list(beta_agg = tmp1, beta = rifle_beta_list, naive_Zs = naive_Zs,  deflated_Zs =  deflated_Zs))
}
