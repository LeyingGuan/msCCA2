# class object for mCCA
my_init = function(xlist, A = 4){
  D = length(xlist)
  ps = sapply(xlist, function(z) dim(z)[2])
  pss = c(0,cumsum(ps))
  ptotal = pss[D+1]
  n = nrow(xlist[[1]])
  xagg = matrix(0, nrow = n, ncol = ptotal)
  for(d in 1:D){
    ll = (pss[d]+1):pss[d+1]
    xagg[,ll] = xlist[[d]]
  }
  S = t(xagg)%*%xagg/n
  Lambda = array(0, dim = dim(S))

  for(d in 1:D){
    ll = (pss[d]+1):pss[d+1]
    Lambda[ll,ll] = t(xagg[,ll])%*%xagg[,ll]/n
    S[ll,ll] = 0
  }
  #find a best threshold m
  ms = ceiling(seq(from = sqrt(n/log(ptotal)), to = n/log(ptotal), length.out = 3))
  rho_tr = rep(0, length(ms))
  beta_inits_all = list()
  beta_inits_sub = list()
  for(i in 1:length(ms)){
    xlist_sub = list()
    S1 = S
    m =ms[i]
    thr = sort(abs(as.vector(S)),decreasing = T)[m^2]
    S1 = sign(S1) * ifelse(abs(S1)<thr,0, abs(S1)-thr)
    ss = apply(S^2,1,sum)
    #take top n/D features from each assay
    n0 = ceiling(n/(  A*D))
    idx_block = list()
    idx_combined = c()
    for(d in 1:D){
      ll = (pss[d]+1):pss[d+1]
      ss0 = ss[ll]
      thr0 = sort(ss0,decreasing = T)[n0]
      idx_block[[d]] = which(ss0>=thr0)
      xlist_sub[[d]] = xlist[[d]][,idx_block[[d]]]
      idx_combined = c(idx_combined, ll[idx_block[[d]]])
    }
    tmp1 =  try(rgcca(A = xlist_sub,  tau = "optimal",
                      scheme = "horst", verbose =F, ncomp = rep(1,length(xlist))))
    if("try-error"%in%class(tmp1)){
      for(d in 1:D){
        beta_inits_sub[[i]][[d]] =  rnorm(ps[d])
      }
    }else{
      beta_inits_sub[[i]] = tmp1$astar
    }
    #achived tau
    beta_inits_all[[i]] = list()
    for(d in 1:D){
      beta_inits_all[[i]][[d]] = rep(0,ps[d])
      beta_inits_all[[i]][[d]][idx_block[[d]]] = tmp1$astar[[d]]
    }
    tmp2 = rho_estimation(xlist_sub, beta_inits_sub[[i]])
    rho_tr[i] = tmp2$rho[1]
    # print("##########")
    # print(cor(Zsum.truth,tmp2$Zsum))
  }
  idx = which.max(rho_tr)
  beta_init0 = beta_inits_all[[idx]]
  return(beta_init0)
}

#' mCCA: msCCA object
#'@import R6
#'@import Rcpp
#'@useDynLib solvers, .registration=TRUE
#'@export
#'\examples{
#'}
mCCA = R6::R6Class(classname = "msCCAobj",public= list(
  norm_type = NULL,
  rho = NULL,
  X = NULL,
  R = NULL,
  beta_init = NULL,
  eta = NULL,
  cl = NULL,
  eta_ratio = NULL,
  eta_low = NULL,
  eps = NULL,
  D = NULL,
  n = NULL,
  ps = NULL,
  pss = NULL,
  ptotal = NULL,
  print_out = NULL,
  rho_tol = NULL,
  rho_maxit = NULL,
  l1proximal_tol = NULL,
  l1proximal_maxit = NULL,
  line_search = NULL,
  line_maxit = NULL,
  out_full = NULL,
  out_cv_tmp = NULL,
  out_cv = NULL,
  norm_grids = NULL,
  Zsum_te = NULL,
  Zs_te = NULL,
  rho_te = NULL,
  aic.approx = NULL,
  prev_size = 0,
  prev_contribution = rep(0,1),
  prev_directions_agg = matrix(0,2,2),
  prev_directions = list(),
  prev_Zsum = matrix(0,2,2),
  prev_Zsum_orth = matrix(0,2,2),
  prev_norms = 0,
  init_method = NULL, 
  trace =NULL,
  my_init_param = NULL,
  initialize = function(X, beta_init = NULL, norm_type = 1, cl = 0.1, eta = 0.2, eta_ratio = NULL,
                        rho_tol = 1e-6, rho_maxit = 1e3, l1proximal_tol =1e-4, l1proximal_maxit = 1e3,
                        line_search = TRUE, line_maxit = 10, eta_low = NULL, eps = NULL,
                        init_method = "rgcca", my_init_param = NULL,
                        trace = TRUE, print_out = 50){
    self$X = X;
    self$R = X;
    self$trace = trace;
    self$norm_type = norm_type;
    self$D = length(X)
    self$n = nrow(X[[1]])
    self$ps = sapply(self$X, function(z) dim(z)[2])
    pss = cumsum(self$ps)
    self$pss = c(0, pss)
    self$ptotal = pss[self$D]
    self$init_method = init_method
    if(is.null(beta_init)){
      if(is.null(my_init_param)){
        self$my_init_param = 4
      }else{
        self$my_init_param = my_init_param
      }
      self$beta_init = self$beta_init_func(self$R)
    }else{
      self$beta_init = beta_init;
    }
    if(is.null(eta_ratio)){
      eta_ratio = 1.0/log(self$n, base = 10)
    }
    if(is.null(eta_low)){
      eta_low = 1.0/self$n
    }
    if(is.null(eps)){
      eps = log(self$ptotal)/self$n
    }
    self$eps = eps
    self$cl = cl
    self$eta = eta
    self$eta_ratio = eta_ratio
    self$rho_tol= rho_tol
    self$rho_maxit = rho_maxit
    self$l1proximal_tol = l1proximal_tol
    self$l1proximal_maxit = l1proximal_maxit
    self$line_search = line_search
    self$line_maxit = line_maxit
    self$eta_low = eta_low
    self$print_out = print_out
  },
  beta_init_func = function(R){
    if(self$init_method == "rgcca"){
      tmp<- try(rgcca(A =R,  tau = "optimal",
                      scheme = "horst", verbose =F, ncomp = rep(1,length(R))))
      if(class(tmp)!="try-error"){
        beta_init0 = tmp$astar
      }else{
        beta_init0 = list()
        for(d in 1:self$D){
          beta_init0[[d]] = rnorm(self$ps[d])
        }
      }
    }else if(self$init_method == "pma"){
      penalties = sqrt(self$n/(log(self$ps)))
      penalties=ifelse(penalties < sqrt(self$ps/4), penalties,sqrt(self$ps/4))
      tmp<- try(MultiCCA(R, type=rep("standard", length(R)),
                         penalty=penalties, ncomponents=1,   trace = F))
      if(class(tmp)!="try-error"){
        beta_init0 = tmp$ws
      }else{
        beta_init0 = list()
        for(d in 1:self$D){
          beta_init0[[d]] = rnorm(self$ps[d])
        }
      }
      }else if(self$init_method == "convex"){
        Ragg =  array(0, dim =c(nrow(R[[1]]),self$ptotal))
        Lambda = array(0, dim = c(self$ptotal, self$ptotal))
        for(d in 1:D){
          ll = (self$pss[d]+1):self$pss[d+1]
          Ragg[,ll] = R[[d]]
          Lambda[ll,ll] = t(self$X[[d]])%*%self$X[[d]]/n
        }
        Sigma = t(Ragg)%*%Ragg/nrow(R[[1]])
        tmp <- initial.convex(A = Sigma, B =Lambda , lambda = sqrt(log(self$D)/self$n), K = 1, nu = 1, epsilon = 0.05, maxiter = 100, trace = self$trace)
        beta_init0 = list()
        u = svd(tmp$Pi)$u[,1]
        for(d in 1:self$D){
          ll = (self$pss[d]+1):self$pss[d+1]
          beta_init0[[d]] = u[ll]
        }
        }else if(self$init_method == "soft-thr"){
          beta_init0 = my_init(self$R, self$my_init_param)
          }else{
      beta_init0 = list()
      for(d in 1:self$D){
        beta_init0[[d]] = rnorm(self$ps[d])
      }
    }
    s = 0.0
    for(k in 1:length(beta_init0)){
      s = s+sum(beta_init0[[k]]^2)
    }
    beta_init = list()
    for(k in 1:length(beta_init0)){
      beta_init[[k]] = beta_init0[[k]]/sqrt(s)
    }
    return(beta_init)
  },
  direction_update_l1_single = function(X = NULL, R = NULL,  beta_init = NULL, l1norm = NULL, l0norm = NULL, trace = F){
    if(is.null(l1norm)){
      stop("no l1 norm provided!")
    }
    if(is.null(X)){
      X = self$X
    }
    if(is.null(R)){
      R = self$R
    }
    beta_init1 = list()
    for(d in 1:self$D){
      if(is.null(beta_init)){
        beta_init1[[d]] = self$beta_init[[d]]
      }else{
        beta_init1[[d]]  = beta_init[[d]]
      }
    }
    out =msCCA_proximal_rank1(beta =beta_init1, X = X, R = R, 
                              rho_tol = self$rho_tol , rho_maxit = self$rho_maxit
                              , eta = self$eta, norm_type = 1,  l0norm = 0, l1norm = l1norm, cl =  self$cl, eta_ratio = self$eta_ratio, 
                              l1proximal_tol =  self$l1proximal_tol, l1proximal_maxit =  self$l1proximal_maxit, line_search =  self$line_search, 
                              line_maxit =  self$line_maxit, eta_low =  self$eta_low, eps = self$eps*l1norm^1.5, trace =  trace, print_out =  self$print_out)
    return(out)
  },
  direction_update_l0_single = function(X = NULL, R = NULL, beta_init = NULL, l1norm = NULL, l0norm = NULL, trace = F){
    if(is.null(l0norm)){
      stop("no l0 norm provided!")
    }
    if(is.null(X)){
      X = self$X
    }
    if(is.null(R)){
      R = self$R
    }
    beta_init1 = list()
    for(d in 1:self$D){
      if(is.null(beta_init)){
        beta_init1[[d]] = self$beta_init[[d]]
      }else{
        beta_init1[[d]]  = beta_init[[d]]
      }
    }
    out =msCCA_proximal_rank1(beta =beta_init1, X = X, R = R, 
                              rho_tol = self$rho_tol , rho_maxit = self$rho_maxit
                              , eta = self$eta, norm_type = 0,  l0norm = l0norm, l1norm = 0, cl =  self$cl, eta_ratio = self$eta_ratio, 
                              l1proximal_tol =  self$l1proximal_tol, l1proximal_maxit =  self$l1proximal_maxit, line_search =  self$line_search, 
                              line_maxit =  self$line_maxit, eta_low =  self$eta_low, eps = self$eps*l0norm^(0.75), trace =  trace, print_out =  self$print_out)
    return(out)
  },
  direction_update_grid = function(X = NULL, R = NULL, beta_init = NULL, norm_grids = NULL, trace = F, update = T, silence = F){
    if(is.null(norm_grids)){
      stop("no norm grids provided!")
    }else{
      #check if ordered
      norm_grids = sort(norm_grids, decreasing = T)
    }
    self$norm_grids = norm_grids
    if(is.null(X)){
      X = self$X
    }
    if(is.null(R)){
      R = self$R
    }
    if(is.null(beta_init)){
      beta_init = list()
      for(d in 1:self$D){
        beta_init[[d]] = self$beta_init[[d]]
      }
    }
    M = length(norm_grids)
    betas = list()
    beta_aggs = array(0, dim = c(self$ptotal, M))
    rhos = rep(0, M)
    l1bounds = rep(0, M)
    n = nrow(X[[1]])
    Zs = array(0, dim = c(n, self$D, M))
    Zsum = array(0, dim = c(n, M))
    etas = rep(0, M)
    for(l in 1:M){
      if(!silence){
        print(l)
      }
      if(self$norm_type == 1){
        out = self$direction_update_l1_single(X = X, R = R, beta_init = beta_init, l1norm = norm_grids[l], trace = trace)
      }else{
        out = self$direction_update_l0_single(X = X, R = R, beta_init = beta_init, l0norm = norm_grids[l], trace = trace)
      }
      l1bounds[l] = out$l1bound
      Zs[,,l] = out$Zs
      Zsum[,l] = out$Zsum
      betas[[l]] = out$beta
      beta_aggs[,l] = out$beta_aug
      rhos[l] = out$rho
      etas[l] = out$eta
    }
    res_out = list(betas = betas, beta_aggs = beta_aggs, Zs = Zs, Zsum = Zsum, rhos = rhos,  l1bounds =  l1bounds, norm_grids = norm_grids)
    if(update){
      self$out_full = res_out
    }else{
      return(res_out)
    }
  },
  direction_cv_grid = function(nfolds = 10, foldid = NULL, seed = 2021, n.core = NULL, trace = F){
    set.seed(seed)
    if(is.null(n.core)){
      n.core = min(nfolds, detectCores())
    }
    if(is.null(foldid)){
      tmp = rep(1:nfolds, each = floor(self$n/nfolds))
      if(length(tmp)<self$n){
        tmp = c(c(1:(self$n-length(tmp))),tmp)
      }
      foldid = sample(1:n, self$n, replace = F)
      foldid = tmp[foldid]
    }else{
      if(length(foldid)!=self$n){
        stop("foldid length != n")
      }else{
        aa = sort(unique(foldid))
        bb = c(1:length(aa))
        cc = sum(aa!=bb)
        if(cc!=0){
          stop("foldid must be from continous integer set from 1 to nfolds!")
        }
      }
      nfolds = length(unique(foldid))
    }
    norm_grids =  self$norm_grids
    beta_cv_lists = list()
    Zs_cv = array(0, dim = c(self$n,self$D ,length(norm_grids)))
    Zsum_cv = array(0, dim = c(self$n,length(norm_grids)))
    rho_cv = rep(0, length(norm_grids))
    for(d in 1:self$D){
      beta_cv_lists[[d]] = array(0, dim = c(self$ps[d], nfolds,length(norm_grids)))
    }
    beta_cv_aggs = array(0, dim = c(self$ptotal,length(norm_grids), nfolds))
    fold_ids = sort(unique(foldid))
    fun0 = function(input0, silence = F){
      nf = input0[1]
      silence = ifelse(input0[2]==1,F, T)
      if(trace){
        print(paste0("fold", nf))       
      }

      Xtmp = list()
      Rtmp = list()
      for(d in 1:self$D){
        Xtmp[[d]] =  self$X[[d]][foldid!=nf,]
        Rtmp[[d]] =  self$R[[d]][foldid!=nf,]
      }
      beta_init = self$beta_init_func(Rtmp)
      out_tmp = self$direction_update_grid(X = Xtmp, R = Rtmp, beta_init = beta_init, norm_grids = norm_grids,trace = F, update = F,  silence =  silence)
      return(out_tmp)
    }
    inputs = list()
    for(i in 1:length(fold_ids)){
      if(fold_ids[i]%%n.core==1 & trace){
        inputs[[i]] = c(fold_ids[i], 1)
      }else{
        inputs[[i]] = c(fold_ids[i], 0)
      }
      
    }
    outputs <-try(mclapply(inputs,fun0, mc.cores =n.core))
    if(class(outputs) == "try-error"){
      print(class(outputs))
      for(i in 1:length(fold_ids)){
        inputs[[i]][2] = 1
      }
      outputs = lapply( inputs, fun0)
    }
    for(nf in 1:nfolds){
      for(d in 1:self$D){
        for(l in 1:length(norm_grids)){
          beta_cv_lists[[d]][,nf,l] =  outputs[[nf]]$betas[[l]][[d]]
          beta_cv_aggs[(self$pss[d]+1):(self$pss[d+1]),l,nf] = outputs[[nf]]$betas[[l]][[d]]
          Zs_cv[foldid == nf,d, l] = 1/sqrt(self$n) * self$R[[d]][foldid == nf,] %*%beta_cv_lists[[d]][,nf,l]
          Zsum_cv[foldid==nf,l] = Zsum_cv[foldid==nf,l]+Zs_cv[foldid == nf,d, l]
        }
      }
    }
    for(l in 1:length(norm_grids)){
      rho_cv[l] = sum(Zsum_cv[,l]^2)/sum(Zs_cv[,,l]^2)
    }
    res_cv = list(beta_cv = beta_cv_lists, Zs_cv = Zs_cv, Zsum_cv = Zsum_cv, rho_cv = rho_cv, foldid = foldid)
    self$out_cv = res_cv
  },
  prediction = function(xlist.te){
    n.te = nrow(xlist.te[[1]])
    Zs.te = array(0, dim = c(n.te, D, length(self$norm_grids)))
    Zsum.te = array(0, dim = c(n.te, length(self$norm_grids)))
    rho.te = rep(0, length(self$norm_grids))
    for(l in 1:length(self$norm_grids)){
      for(d in 1:D){
        Zs.te[,d,l] = 1/sqrt(n.te) *xlist.te[[d]]%*%self$out_full$betas[[l]][[d]]
        Zsum.te[,l] = Zsum.te[,l]+Zs.te[,d,l]
      }
    }
    A= apply(Zsum.te^2,2,sum)
    B = apply(Zs.te^2,3,sum)
    self$Zsum_te = Zsum.te
    self$Zs_te = Zs.te
    self$rho_te  = A/B
  },
  plot_cv_func = function(){
    norm_grids = self$norm_grids
    if(!is.null(self$rho_te)){
      ymin =min(c( self$rho_te, self$out_cv$rho_cv, self$out_full$rhos)-0.5)
      ymax = max(c( self$rho_te, self$out_cv$rho_cv, self$out_full$rhos)+0.5)
    }else{
      ymin = min(c(self$out_cv$rho_cv, self$out_full$rhos))
      ymax = max(c(self$out_cv$rho_cv, self$out_full$rhos))
    }
    
    plot(norm_grids, self$out_cv$rho_cv, ylim= c(ymin, ymax), xlab = "norms",ylab = "rho")
    points(norm_grids,msCCAproximal_l1$out_full$rhos, col = "blue", pch = 19)
    if(!is.null(self$rho_te)){
      points(norm_grids, self$rho_te, col = "red", pch = 19)
      legend("bottomright", legend = c("tr","te","cv"),pch = c( 19, 19,1), col = c("blue", "red", "black"),bty = "n")
    }else{
      legend("bottomright", legend = c("tr","cv"),pch = c( 19,1), col = c("blue", "black"),bty = "n")
    }
    abline(v = norm_grids[which.max(self$out_cv$rho_cv)], lty = 2)
  },
  direction_grow = function(model_idx=1){
    self$prev_size = self$prev_size+1
    if(self$prev_size==1){
      self$prev_norms = self$out_full$norm_grids[model_idx]
      self$prev_contribution = self$out_full$rhos[model_idx]
      self$prev_directions_agg  = matrix(self$out_full$beta_aggs[,model_idx],ncol = 1)
      self$prev_directions = list()
      for(d in 1:self$D){
        ii = (self$pss[d]+1):self$pss[d+1]
        tmp = self$X[[d]]%*%self$out_full$betas[[model_idx]][[d]]
        tmp =  t(self$X[[d]])%*%tmp/self$n
        self$prev_directions[[d]] = self$out_full$betas[[model_idx]][[d]]
      }
    }else{
      self$prev_contribution = c(self$prev_contribution, self$out_full$rhos[model_idx])
      self$prev_directions_agg = cbind(self$prev_directions_agg, self$out_full$beta_aggs[,model_idx])
      self$prev_norms = c(self$prev_norms , self$out_full$norm_grids[model_idx])
      tmp =   rep(0, self$ptotal)
      for(d in 1:self$D){
        ii = (self$pss[d]+1):self$pss[d+1]
        tmp1 = self$X[[d]]%*%self$out_full$betas[[model_idx]][[d]]
        tmp[ii] = t(self$X[[d]])%*%tmp1/self$n
      }
      for(d in 1:self$D){
        self$prev_directions[[d]] = cbind(self$prev_directions[[d]],self$out_full$betas[[model_idx]][[d]])
      }
    }
    #update residual for better initialization
    zsum = self$out_full$Zsum[,model_idx]
    zsum_orth  = zsum
    if(self$prev_size == 1){
      self$prev_Zsum = matrix(zsum, ncol = 1)
      self$prev_Zsum_orth = matrix(zsum_orth, ncol = 1)
    }else{
      self$prev_Zsum = cbind(self$prev_Zsum, zsum)
      for(k in 1:(self$prev_size-1)){
        zsum_orth = zsum_orth - self$prev_Zsum_orth[,k] * sum(self$prev_Zsum_orth[,k]* zsum_orth)/sum(self$prev_Zsum_orth[,k]*self$prev_Zsum_orth[,k])
      }
      self$prev_Zsum_orth = cbind(self$prev_Zsum_orth, zsum_orth)
    }
    for(d in 1:self$D){
      tmp =  t(self$R[[d]])%*%zsum_orth/sum(zsum_orth^2)
      self$R[[d]] = self$R[[d]] -zsum_orth%*%t(tmp)
    }
  },
  direction_grow_lazy = function(out, norm){
    self$prev_size = self$prev_size+1
    if(self$prev_size==1){
      self$prev_norms = norm
      self$prev_contribution = out$rho
      self$prev_directions_agg  = matrix(out$beta_aug,ncol = 1)
      self$prev_directions = list()
      for(d in 1:self$D){
        ii = (self$pss[d]+1):self$pss[d+1]
        tmp = self$X[[d]]%*%out$beta[[d]]
        tmp =  t(self$X[[d]])%*%tmp/self$n
        self$prev_directions[[d]] = out$beta[[d]]
      }
    }else{
      self$prev_contribution = c(self$prev_contribution, out$rho)
      self$prev_directions_agg = cbind(self$prev_directions_agg,out$beta_aug)
      self$prev_norms = c(self$prev_norms , norm)
      tmp =   rep(0, self$ptotal)
      for(d in 1:self$D){
        ii = (self$pss[d]+1):self$pss[d+1]
        tmp1 = self$X[[d]]%*%out$beta[[d]]
        tmp[ii] = t(self$X[[d]])%*%tmp1/self$n
      }
      for(d in 1:self$D){
        self$prev_directions[[d]] = cbind(self$prev_directions[[d]],out$beta[[d]])
      }
    }
    #update residual for better initialization
    zsum = out$Zsum
    zsum_orth  = zsum
    if(self$prev_size == 1){
      self$prev_Zsum = matrix(zsum, ncol = 1)
      self$prev_Zsum_orth = matrix(zsum_orth, ncol = 1)
    }else{
      self$prev_Zsum = cbind(self$prev_Zsum, zsum)
      for(k in 1:(self$prev_size-1)){
        zsum_orth = zsum_orth - self$prev_Zsum_orth[,k] * sum(self$prev_Zsum_orth[,k]* zsum_orth)/sum(self$prev_Zsum_orth[,k]*self$prev_Zsum_orth[,k])
      }
      self$prev_Zsum_orth = cbind(self$prev_Zsum_orth, zsum_orth)
    }
    for(d in 1:self$D){
      tmp =  t(self$R[[d]])%*%zsum_orth/sum(zsum_orth^2)
      self$R[[d]] = self$R[[d]] -zsum_orth%*%t(tmp)
    }
  }
)
)


#' msCCAl1: wrapper function for msCCA with decaying l1 norm bound for multiple components
#'@export
#'\examples{
#'}
msCCAl1 = function(xlist, ncomp, xlist.te =NULL, init_method = "soft-thr", nfolds = 7, foldid = NULL,
                   norms = NULL, norm_grids = NULL, nlambda = 10,
                   rho_maxit = 200, print_out = 100, l1proximal_maxit = 1e4, rho_tol = 1e-2,
                   eps =NULL, mode = "full",trace = F){
  D = length(xlist)
  ps = sapply(xlist,function(z) dim(z)[2])
  pss = c(0,cumsum(ps))
  ptotal = pss[D+1]
  n = nrow(xlist[[1]])
  if(is.null(norms)){
    if(is.null(norm_grids)){
      norm_grids = sort(exp(seq(from = log(sqrt(5.0)), to =log(min(sqrt(ptotal),sqrt(n/2))), length.out =nlambda)), decreasing = T)
    }
  }else{
    if(length(norms)!=ncomp){
      stop("length of norms != ncomp.")
    }
  }
  if(is.null(foldid)){
    foldid = sample(rep(1:nfolds, each = ceiling(n/nfolds)), n)
  }else{
    nfolds = length(unique(foldid))
  }
  
  if(is.null(eps)){
    eps = log(ptotal)/n
  }
  msCCAproximal_l1 = mCCA$new(X = xlist, beta_init =NULL, norm_type = 1, init_method = init_method,
                              rho_maxit = rho_maxit, print_out = print_out,
                              l1proximal_maxit = l1proximal_maxit, rho_tol = rho_tol,
                              eps = eps, trace = trace)
  beta_inits = list()
  selected_penalties = rep(0, ncomp)
  for(k in 1:ncomp){
    print(paste0("####################comp", k))
    beta_inits[[k]] = msCCAproximal_l1$beta_init
    if(!is.null(norms)){
      out = msCCAproximal_l1$direction_update_l1_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R, 
                                                        beta_init =  msCCAproximal_l1$beta_init, l1norm = norms[k], trace = trace)
      #update
      selected_penalties[k] = norms[k]
      msCCAproximal_l1$direction_grow_lazy(out, norm = norms[k])
    }else{
      if(mode == "lazy"){
        #create a norm with equal value
        msCCAproximal_l1$direction_update_grid(norm_grids = norm_grids, trace = trace , update = T)
        msCCAproximal_l1$direction_cv_grid(foldid = foldid, trace = trace)
        msCCAproximal_l1$out_cv$rho_cv[is.na(msCCAproximal_l1$out_cv$rho_cv)] = 0
        selected_penalties[k] =  which.max(msCCAproximal_l1$out_cv$rho_cv)
        msCCAproximal_l1$direction_grow(model_idx =selected_penalties[k])
        selected_penalties[k] = norm_grids[selected_penalties[k]]
        norms = rep(selected_penalties[k], ncomp)
      }else{
        #do not create norms and tune every time
        msCCAproximal_l1$direction_update_grid(norm_grids = norm_grids, trace = trace , update = T)
        msCCAproximal_l1$direction_cv_grid(foldid = foldid, trace = trace)
        msCCAproximal_l1$out_cv$rho_cv[is.na(msCCAproximal_l1$out_cv$rho_cv)] = 0
        selected_penalties[k] =  which.max(msCCAproximal_l1$out_cv$rho_cv)
        msCCAproximal_l1$direction_grow(model_idx =selected_penalties[k])
        selected_penalties[k] = norm_grids[selected_penalties[k]]
      }
    }
    if(k < ncomp){
      msCCAproximal_l1$beta_init = msCCAproximal_l1$beta_init_func(msCCAproximal_l1$R)
    }
    
  }
  #naive
  aa = rho_estimation(xlist,msCCAproximal_l1$prev_directions); 
  naive_rho.tr = aa$rho; 
  Zs = aa$Z; Zsum = aa$Zsum;
  naive_rho.te = NULL;Zs.te = NULL; Zsum.te = NULL
  if(!is.null(xlist.te)){
    aa = rho_estimation(xlist.te,msCCAproximal_l1$prev_directions); 
    naive_rho.te = aa$rho; 
    Zs.te = aa$Z; Zsum.te = aa$Zsum;
  }
  aa = rho_estimation_deflated(xlist,msCCAproximal_l1$prev_directions); 
  deflated_rho.tr = aa$rho; 
  Zsum.deflated =aa$Zsum;
  Zs.deflated = aa$Z
  deflated_rho.te = NULL
  Zsum.te.deflated = NULL
  Zs.te.deflated = NULL
  if(!is.null(xlist.te)){
    aa = rho_estimation_deflated(xlist.te,msCCAproximal_l1$prev_directions); 
    deflated_rho.te = aa$rho; 
    Zs.te.deflated = aa$Z; Zsum.te.deflated = aa$Zsum;
  }
  naive_Zs = list(Zs = Zs, Zsum = Zsum, Zs.te = Zs.te, Zsum.te = Zsum.te,
                  rho.tr = naive_rho.tr, rho.te = naive_rho.te)
  deflated_Zs = list(Zs = Zs.deflated, Zsum = Zsum.deflated, Zs.te = Zs.te.deflated, 
                     Zsum.te = Zsum.te.deflated,
                     rho.tr = deflated_rho.tr, rho.te = deflated_rho.te)
  #calculate Zsum
  return(list(fitted_model = msCCAproximal_l1, beta_inits = beta_inits, selected_penalties = selected_penalties,  
              naive_Zs = naive_Zs,  deflated_Zs =  deflated_Zs))
}

#' msCCAl0: wrapper function for msCCA with decaying l0 norm bound for multiple components
#'@export
#'\examples{
#'}
msCCAl0 = function(xlist, ncomp, xlist.te =NULL, init_method = "soft-thr", nfolds = 7, foldid = NULL, 
                   norms = NULL, norm_grids = NULL, nlambda = 10,
                   rho_maxit = 200, print_out = 100, rho_tol = 1e-2,
                   eps =NULL, mode = "full",trace = F){
  D = length(xlist)
  ps = sapply(xlist,function(z) dim(z)[2])
  pss = c(0,cumsum(ps))
  ptotal = pss[D+1]
  n = nrow(xlist[[1]])
  if(is.null(norms)){
    if(is.null(norm_grids)){
      norm_grids = floor(sort(exp(seq(from =log(5.0), to =log(min(ptotal,n/2)), length.out =nlambda)), decreasing = T))
    }
  }else{
    if(length(norms)!=ncomp){
      stop("length of norms != ncomp.")
    }
  }
  if(is.null(foldid)){
    foldid = sample(rep(1:nfolds, each = ceiling(n/nfolds)), n)
  }else{
    nfolds = length(unique(foldid))
  }
  if(is.null(eps)){
    eps = log(ptotal)/n
  }
  msCCAproximal_l0 = mCCA$new(X = xlist, beta_init =NULL, norm_type = 0, init_method = init_method,
                              rho_maxit = rho_maxit, print_out = print_out, rho_tol = rho_tol,
                              eps = eps, trace = trace)
  beta_inits = list()
  selected_penalties = rep(0, ncomp)
  for(k in 1:ncomp){
    print(paste0("####################comp", k))
    beta_inits[[k]] = msCCAproximal_l0$beta_init
    
    if(!is.null(norms)){
      out = msCCAproximal_l0$direction_update_l0_single(X =  msCCAproximal_l0$X, R =  msCCAproximal_l0$R, 
                                                        beta_init =  msCCAproximal_l0$beta_init, l0norm = norms[k], trace = trace)
      #update
      selected_penalties[k] = norms[k]
      msCCAproximal_l0$direction_grow_lazy(out, norm = norms[k])
    }else{
      if(mode == "lazy"){
        #create a norm with equal value
        msCCAproximal_l0$direction_update_grid(norm_grids = norm_grids, trace = trace , update = T)
        msCCAproximal_l0$direction_cv_grid(foldid = foldid, trace = trace)
        msCCAproximal_l0$out_cv$rho_cv[is.na(msCCAproximal_l0$out_cv$rho_cv)] = 0
        selected_penalties[k] =  which.max(msCCAproximal_l0$out_cv$rho_cv)
        msCCAproximal_l0$direction_grow(model_idx =selected_penalties[k])
        selected_penalties[k] = norm_grids[selected_penalties[k]]
        norms = rep(selected_penalties[k], ncomp)
      }else{
        #do not create norms and tune every time
        msCCAproximal_l0$direction_update_grid(norm_grids = norm_grids, trace = trace , update = T)
        msCCAproximal_l0$direction_cv_grid(foldid = foldid, trace = trace)
        msCCAproximal_l0$out_cv$rho_cv[is.na(msCCAproximal_l0$out_cv$rho_cv)] = 0
        selected_penalties[k] =  which.max(msCCAproximal_l0$out_cv$rho_cv)
        msCCAproximal_l0$direction_grow(model_idx =selected_penalties[k])
        selected_penalties[k] = norm_grids[selected_penalties[k]]
      }
    }
    if(k < ncomp){
      msCCAproximal_l0$beta_init = msCCAproximal_l0$beta_init_func(msCCAproximal_l0$R)
    }
    
  }
  #naive
  aa = rho_estimation(xlist,msCCAproximal_l0$prev_directions); 
  naive_rho.tr = aa$rho; 
  Zs = aa$Z; Zsum = aa$Zsum;
  naive_rho.te = NULL;Zs.te = NULL; Zsum.te = NULL
  if(!is.null(xlist.te)){
    aa = rho_estimation(xlist.te,msCCAproximal_l0$prev_directions); 
    naive_rho.te = aa$rho; 
    Zs.te = aa$Z; Zsum.te = aa$Zsum;
  }
  aa = rho_estimation_deflated(xlist,msCCAproximal_l0$prev_directions); 
  deflated_rho.tr = aa$rho; 
  Zsum.deflated =aa$Zsum;
  Zs.deflated = aa$Z
  deflated_rho.te = NULL
  Zsum.te.deflated = NULL
  Zs.te.deflated = NULL
  if(!is.null(xlist.te)){
    aa = rho_estimation_deflated(xlist.te,msCCAproximal_l0$prev_directions); 
    deflated_rho.te = aa$rho; 
    Zs.te.deflated = aa$Z; Zsum.te.deflated = aa$Zsum;
  }
  naive_Zs = list(Zs = Zs, Zsum = Zsum, Zs.te = Zs.te, Zsum.te = Zsum.te,
                  rho.tr = naive_rho.tr, rho.te = naive_rho.te)
  deflated_Zs = list(Zs = Zs.deflated, Zsum = Zsum.deflated, Zs.te = Zs.te.deflated, 
                     Zsum.te = Zsum.te.deflated,
                     rho.tr = deflated_rho.tr, rho.te = deflated_rho.te)
  #calculate Zsum
  return(list(fitted_model = msCCAproximal_l0, beta_inits = beta_inits, selected_penalties = selected_penalties,  
              naive_Zs = naive_Zs,  deflated_Zs =  deflated_Zs))
}





