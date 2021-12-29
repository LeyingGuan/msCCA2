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
      thr0 = sort(ss0,decreasing = T)[max(min(n0,length(ss0)),1)]
      idx_block[[d]] = which(ss0>=thr0)
      xlist_sub[[d]] = xlist[[d]][,idx_block[[d]]]
      idx_combined = c(idx_combined, ll[idx_block[[d]]])
    }
    tmp1 =  try(rgcca(A = xlist_sub,  tau = "optimal",
                      scheme = "horst", verbose =F, ncomp = rep(1,length(xlist_sub))))
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


deflation_func = function(Ragg, R, zsum_orth){
  Ragg = lm(Ragg~zsum_orth)$residuals
  D = length(R)
  ps = sapply(R, function(z) dim(z)[2])
  pss = cumsum(ps)
  pss = c(0, pss)
  Rnew = list()
  for(d in 1:D){
    ll = (pss[d]+1):pss[d+1]
    Rnew[[d]] = Ragg[,ll]
  }
  return(list(Ragg = Ragg, R = Rnew))
}


zsum_func = function(Xagg, beta_agg, prev_Zsum_orth = NULL){
  #update residual for better initialization
  n = nrow(Xagg)
  zsum = Xagg%*%beta_agg/sqrt(n)
  zsum_orth = zsum
  if(!is.null(prev_Zsum_orth)){
    zsum_orth = lm(zsum~prev_Zsum_orth)$residuals
  }
  return(list(zsum = zsum, zsum_orth = zsum_orth))
}

direction_single_update_l1_func = function(X, R,  beta_init, eta, eta_ratio,eta_low,  eps, 
                                        l1norm_max , l1norm_min,  rho_tol, rho_maxit,
                                        l1proximal_tol, l1proximal_maxit, line_maxit, 
                                        warm_up, trace, print_out, early_stop, norm_varying_ridge){
  D = length(X)
  out =msCCA_proximal_rank1(beta =beta_init, X = X, R = R, eta = eta, eta_ratio = eta_ratio,
                            eta_low =  eta_low,   eps = eps, l1norm_max = l1norm_max, l1norm_min = l1norm_min,
                            rho_tol =  rho_tol , rho_maxit = rho_maxit,  norm_type = 1,  l0norm = 0, 
                             l1proximal_tol = l1proximal_tol,  l1proximal_maxit =  l1proximal_maxit,
                            line_maxit =  line_maxit,warm_up = warm_up, trace =  trace, print_out =  print_out,
                            early_stop = early_stop, norm_varying_ridge=norm_varying_ridge)
  names = c("beta_augs", "beta_norms", "bounds",  "rho_denominators", "rho_numerators", "rhos")
  for(i in 1:length(names )){
    if(names[i]=="beta_augs"){
      out[[names[i]]] = out[[names[i]]][,1:(out$end_iter+1)]
    }else{
      out[[names[i]]] = out[[names[i]]][1:(out$end_iter+1)]
    }
    
  }
  return(out)
  
}


beta_init_func = function(R, init_method, my_init_param=NULL, X = NULL){
  D = length(R)
  n = nrow(R[[1]])
  ps = sapply(R, function(z) ncol(z))
  ptotal = sum(ps)
  pss = c(0, cumsum(ps))
  if(init_method == "rgcca"){
    tmp<- try(rgcca(A =R,  tau = "optimal",
                    scheme = "horst", verbose =F, ncomp = rep(1,length(R))))
    if(class(tmp)!="try-error"){
      beta_init0 = tmp$astar
    }else{
      beta_init0 = list()
      for(d in 1:D){
        beta_init0[[d]] = rnorm(ps[d])
      }
    }
  }else if(init_method == "pma"){
    penalties = sqrt(n/(log(ps)))
    penalties=ifelse(penalties < sqrt(ps/4), penalties,sqrt(ps/4))
    tmp<- try(MultiCCA(R, type=rep("standard", length(R)),
                       penalty=penalties, ncomponents=1,   trace = F))
    if(class(tmp)!="try-error"){
      beta_init0 = tmp$ws
    }else{
      beta_init0 = list()
      for(d in 1:D){
        beta_init0[[d]] = rnorm(ps[d])
      }
    }
  }else if(init_method == "convex"){
    Ragg =  array(0, dim =c(nrow(R[[1]]),ptotal))
    Lambda = array(0, dim = c(ptotal, ptotal))
    for(d in 1:D){
      ll = (pss[d]+1):pss[d+1]
      Ragg[,ll] = R[[d]]
      Lambda[ll,ll] = t(X[[d]])%*%X[[d]]/n
    }
    Sigma = t(Ragg)%*%Ragg/nrow(R[[1]])
    tmp <- initial.convex(A = Sigma, B =Lambda , lambda = sqrt(log(D)/n), K = 1, nu = 1, epsilon = 0.05, maxiter = 1000, trace = F)
    beta_init0 = list()
    u = svd(tmp$Pi)$u[,1]
    for(d in 1:D){
      ll = (pss[d]+1):pss[d+1]
      beta_init0[[d]] = u[ll]
    }
  }else if(init_method == "soft-thr"){
    beta_init0 = my_init(R, my_init_param)
  }else{
    beta_init0 = list()
    for(d in 1:D){
      beta_init0[[d]] = rnorm(ps[d])
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
}
#' msCCAl1: R6 msCCAl1 object: users can use this object to estimate mCCA direction via msCCAl1 sequentially,
#'  and can grow new directions as needed. This is the most flexible way of using msCCAl1 for experienced users.
#'@export
#'\examples{
#'}
#'
#'
msCCAl1 = R6::R6Class(classname = "msCCAl1obj",public= list(
  norm_type = NULL,
  rho = NULL,
  X = NULL,
  R = NULL,
  Xagg = NULL,
  Ragg = NULL,
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
  Zsum_te = NULL,
  Zs_te = NULL,
  rho_te = NULL,
  penalty.C = NULL,
  cv_rhos = NULL,
  penalized_rhos = NULL,
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
  out_single_update = NULL,
  warm_up = NULL,
  l1norm_max = NULL,
  l1norm_min = NULL,
  norm_varying_ridge = NULL,
  initialize = function(X, beta_init = NULL, norm_type = 1, eta = NULL, eta_ratio = NULL,
                        rho_tol = 1e-3, rho_maxit = 1e3, l1proximal_tol =1e-4, l1proximal_maxit = 1e3,
                        line_search = TRUE, line_maxit = 10, eta_low = NULL, eps = NULL,
                        init_method = "rgcca", my_init_param = NULL, norm_varying_ridge = T,
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
    self$norm_varying_ridge = norm_varying_ridge
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
      eta_ratio = sqrt(1/self$n)
    }
    if(is.null(eta_low)){
      eta_low = 1/self$n
    }
    if(is.null(eps)){
      eps = log(self$ptotal)/self$n
    }
    if(is.null(eta)){
      eta = sqrt(1/self$n)
    }
    self$eps = eps
    self$eta = eta
    self$eta_ratio = eta_ratio
    self$rho_tol= rho_tol
    self$rho_maxit = rho_maxit
    self$l1proximal_tol = l1proximal_tol
    self$l1proximal_maxit = l1proximal_maxit
    self$line_maxit = line_maxit
    self$eta_low = eta_low
    self$print_out = print_out
    self$Xagg = self$X[[1]]
    self$Ragg = self$R[[1]]
    for(d in 2:self$D){
      self$Xagg = cbind(self$Xagg, self$X[[d]])
      self$Ragg = cbind(self$Ragg, self$R[[d]])
    }
  },
  beta_init_func = function(R){
    beta_init = beta_init_func(R = R, init_method = self$init_method, my_init_param=self$my_init_param, X = self$X)
    return(beta_init)
  },
  direction_update_single = function(X = NULL, R = NULL,  beta_init = NULL, l1norm_max = NULL, 
                                        l1norm_min = NULL,  warm_up = 50, rho_maxit = NULL, trace = F){
    if((is.null(l1norm_max) | is.null(l1norm_min))){
      stop("no l1 norm max or min provided!")
    }
    self$l1norm_max = l1norm_max
    self$l1norm_min = l1norm_min
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
    self$warm_up = warm_up
    if(trace){
      print(self$warm_up)
    }
    out = direction_single_update_l1_func(X = X, R = R,  beta_init = beta_init1, eta = self$eta, eta_ratio =  self$eta_ratio,
                                          eta_low = self$eta_low,  eps = self$eps,  l1norm_max = l1norm_max , l1norm_min = l1norm_min,  
                                          rho_tol = self$rho_tol, rho_maxit = self$rho_maxit,
                                           l1proximal_tol = self$l1proximal_tol, l1proximal_maxit =  self$l1proximal_maxit, line_maxit =self$line_maxit, 
                                           warm_up = self$warm_up, trace = trace, print_out = self$print_out, early_stop = T,
                                           norm_varying_ridge = self$norm_varying_ridge)
    #   
    # out =msCCA_proximal_rank1(beta =beta_init1, X = X, R = R, rho_tol = self$rho_tol , rho_maxit = rho_maxit, 
    #                             eta = self$eta, norm_type = 1,  l0norm = 0, l1norm_max = l1norm_max, l1norm_min = l1norm_min,
    #                             eta_ratio = self$eta_ratio, l1proximal_tol =  self$l1proximal_tol, 
    #                             l1proximal_maxit =  self$l1proximal_maxit,
    #                             line_maxit =  self$line_maxit, eta_low =  self$eta_low, eps = self$eps,warm_up = warm_up, trace =  trace, print_out =  self$print_out,
    #                             early_stop = early_stop)
    names = c("beta_augs", "beta_norms", "bounds",  "rho_denominators", "rho_numerators", "rhos")
    for(i in 1:length(names )){
      if(names[i]=="beta_augs"){
        out[[names[i]]] = out[[names[i]]][,1:(out$end_iter+1)]
      }else{
        out[[names[i]]] = out[[names[i]]][1:(out$end_iter+1)]
      }
      
    }
    self$out_single_update = out
    return(out)
  },
  
  direction_selection = function(method = "penalized objective", penalty.C = 2,nfolds = 10, foldid = NULL, seed = 2021, multi.core = T, n.core = NULL, trace = F){
    if(method == "penalized objective"){
      Zs = array(0, dim = c(n, self$D, dim(self$out_single_update$beta_augs)[2]))
      Zs_residual = array(0, dim = c(n, self$D, dim(self$out_single_update$beta_augs)[2]))
      for(d in 1:self$D){
        ll = c((self$pss[d]+1):self$pss[d+1])
        Zs[,d,] = self$Xagg[,ll]%*%self$out_single_update$beta_augs[ll,]
        Zs_residual[,d,] = self$Ragg[,ll]%*%self$out_single_update$beta_augs[ll,]
      }
      Zsums = apply(Zs_residual, c(1,3), sum)
      rho.hat = apply(Zsums^2,c(2),sum)/apply(Zs^2,c(3),sum)
      evaluation_obj = rho.hat - penalty.C*sqrt(log(self$ptotal)/n)*self$out_single_update$beta_norms
      # idx_tmp= which.max(evaluation_obj)
      # rho_tmp = rho.hat[idx_tmp]
      # evaluation_obj = rho.hat - (1+(penalty.C-1)*((rho_tmp-1)/(self$D-1)))*sqrt(log(ptotal)/n)*self$out_single_update$beta_norms
      outputs = NULL
      self$penalized_rhos = evaluation_obj
      self$penalty.C = penalty.C
    }else{
      #cross validation
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
      X = self$X;R = self$R; pss = self$pss
      eta = self$eta;eta_ratio = self$eta_ratio; eta_low = self$eta_low; eps = self$eps;
      l1norm_max = self$l1norm_max ; l1norm_min = self$l1norm_min; rho_tol = self$rho_tol; rho_maxit = ncol(self$out_single_update[["beta_augs"]]);
      l1proximal_tol = self$l1proximal_tol;l1proximal_maxit =  self$l1proximal_maxit;line_maxit =self$line_maxit;
      warm_up = self$warm_up;  print_out = self$print_out; init_method = self$init_method; my_init_param = self$my_init_param; norm_varying_ridge =self$norm_varying_ridge
      cv_evaluation = function(fold_id){
        ###use the grid of bounds and step sizes from the full fit
        print(paste0("start fold id = ", fold_id))
        test_foldid = which(foldid==fold_id)
        train_foldid = which(foldid!=fold_id)
        Xtmp = list()
        Rtmp = list()
        D = length(X)
        for(d in 1:D){
          Xtmp[[d]] = X[[d]][train_foldid,]
          Rtmp[[d]] = R[[d]][train_foldid,] 
        }
        beta_init_cv=beta_init_func(R =Rtmp, init_method = init_method , my_init_param=my_init_param, X = Xtmp)
        print(paste0("finish initializing fold id = ", fold_id))
        out_cv = direction_single_update_l1_func(X = Xtmp, R = Rtmp,  beta_init = beta_init_cv, eta = eta, eta_ratio =  eta_ratio,
                                                 eta_low = eta_low,  eps = eps,  l1norm_max = l1norm_max , l1norm_min = l1norm_min,  
                                                 rho_tol = rho_tol, rho_maxit = rho_maxit,l1proximal_tol = l1proximal_tol, 
                                                 l1proximal_maxit =  l1proximal_maxit, line_maxit =line_maxit, 
                                                 warm_up = warm_up, trace = F, print_out = print_out, early_stop = F, norm_varying_ridge = norm_varying_ridge)
        ##evaluate on test data
        Zs = array(NA, dim = c(length(test_foldid), D, dim(out_cv$beta_augs)[2]))
        Zs_residual = array(NA, dim = c(length(test_foldid), D, dim(out_cv$beta_augs)[2]))
        for(d in 1:D){
          ll = c((pss[d]+1):pss[d+1])
          Zs[,d,] =X[[d]][test_foldid,]%*%out_cv$beta_augs[ll,]
          Zs_residual[,d,] =R[[d]][test_foldid,]%*%out_cv$beta_augs[ll,]
        }
        Zsums = apply(Zs_residual, c(1,3), sum)
        B = apply(Zs^2, 3, sum)
        A = apply(Zsums^2, 2, sum)
        rho_components = data.frame(matrix(0, nrow = length(A), ncol = 2))
        colnames(rho_components) = c("numerator", "denominator")
        rho_components[,1] = A
        rho_components[,2] = B
        print(paste0("end fold id = ", fold_id))
        return(rho_components)
      }
      if(multi.core=="mclapply"){
        outputs <-mclapply(1:nfolds,cv_evaluation, mc.cores =n.core)
      }else if(multi.core =="doparallel"){
        #cl <- makePSOCKcluster(n.core)
        #registerDoParallel(n.core)
        outputs <-foreach (i=1:nfolds)%dopar%{
          cv_evaluation(i)
        }

      }else{
        outputs <-lapply(1:nfolds,cv_evaluation)
      }
      print(paste0("finish evaluation "))
      evaluation_obj = outputs[[1]]
      for(i in 2:length(outputs)){
        evaluation_obj= evaluation_obj+outputs[[i]]
      }
      evaluation_obj = evaluation_obj[,1]/evaluation_obj[,2]
      self$cv_rhos = evaluation_obj
    }
    return(list(evaluation_obj = evaluation_obj, outputs = outputs))
  },
  direction_grow = function(step_idx=1){
    self$prev_size = self$prev_size+1
    if(self$prev_size==1){
      self$prev_norms = self$out_single_update$beta_norms[step_idx]
      self$prev_contribution = self$out_single_update$rho_numerators[step_idx]/self$out_single_update$rho_denominators[step_idx]
      self$prev_directions_agg  = matrix(self$out_single_update$beta_augs[,step_idx],ncol = 1)
      self$prev_directions = list()
      for(d in 1:self$D){
        ii = (self$pss[d]+1):self$pss[d+1]
        self$prev_directions[[d]] = self$out_single_update$beta_augs[ii,step_idx]
      }
    }else{
      self$prev_contribution = c(self$prev_contribution, self$out_single_update$rho_numerators[step_idx]/self$out_single_update$rho_denominators[step_idx])
      self$prev_directions_agg = cbind(self$prev_directions_agg,self$out_single_update$beta_augs[,step_idx])
      self$prev_norms = c(self$prev_norms, self$out_single_update$beta_norms[step_idx])
      for(d in 1:self$D){
        ii = (self$pss[d]+1):self$pss[d+1]
        self$prev_directions[[d]] = cbind(self$prev_directions[[d]],self$out_single_update$beta_augs[ii,step_idx])
      }
    }
    #update residual for better initialization
    if(self$prev_size == 1){
      zsum_out = zsum_func(Xagg = self$Xagg, beta_agg = self$prev_directions_agg, prev_Zsum_orth = NULL)
      self$prev_Zsum = matrix(zsum_out$zsum, ncol = 1)
      self$prev_Zsum_orth = matrix(zsum_out$zsum_orth, ncol = 1)
    }else{
      zsum_out = zsum_func(Xagg = self$Xagg, beta_agg = self$prev_directions_agg[,self$prev_size], prev_Zsum_orth = self$prev_Zsum_orth)
      self$prev_Zsum = cbind(self$prev_Zsum, zsum_out$zsum)
      self$prev_Zsum_orth = cbind(self$prev_Zsum_orth, zsum_out$zsum_orth)
    }
    deflate_out = deflation_func(Ragg = self$Ragg, R = self$R, zsum_orth = zsum_out$zsum_orth)
    self$R = deflate_out$R
    self$Ragg = deflate_out$Ragg

  }
  
)
)


#' msCCAl1func: wrapper function for using msCCAl1 to estimate mCCA directions with a pre-determined component number
#' 
#'@export
#'\examples{
#'}
#'
#'

  
msCCAl1func = function(xlist, ncomp, xlist.te = NULL, init_method = "soft-thr", nfolds = 20, warm_up = 50, penalty.C=2, foldid = NULL,
                       l1norm_max = NULL, l1norm_min = NULL,  eta_ratio = NULL, eta = NULL, eps = NULL, my_init_param = NULL, 
                       l1proximal_maxit = 1e4, rho_tol = 1e-3, rho_maxit = 5000, print_out = 100, 
                       step_selection = "penalized", norm_varying_ridge = T,
                       multi.core = "mclapply", seed = 2021, trace = T){
  set.seed(seed)
  D = length(xlist)
  ps = sapply(xlist,function(z) dim(z)[2])
  pss = c(0,cumsum(ps))
  ptotal = pss[D+1]
  n = nrow(xlist[[1]])
  if(is.null(l1norm_max)){
    l1norm_max = min(sqrt(n/4),sqrt(ptotal))
  }
  if(is.null(l1norm_min)){
    l1norm_min = sqrt(2)
  }
  if(is.null(foldid)){
    foldid = sample(rep(1:nfolds, each = ceiling(n/nfolds)), n)
  }else{
    nfolds = length(unique(foldid))
  }
  
  msCCAproximal_l1 = msCCAl1$new(X = xlist, beta_init =NULL, norm_type = 1, init_method = init_method, my_init_param = my_init_param,
                                 eta = eta, eta_ratio = eta_ratio,eps = eps, 
                                 l1proximal_maxit = l1proximal_maxit, rho_tol = rho_tol, rho_maxit = rho_maxit,  norm_varying_ridge=norm_varying_ridge,
                                 print_out = print_out, 
                                 trace = trace)
  
  step_idxs = rep(NA, ncomp)
  l1bounds = rep(NA, ncomp)
  errors_track = list()
  Zsum_orth.te = NA
  if(!is.null(xlist.te)){
    Zsum_orth.te = matrix(0, nrow= nrow(xlist.te[[1]]), ncol = ncomp)
  }
  
  beta_inits = list()
  errors_track_selected = data.frame(matrix(NA, nrow = ncomp, ncol = 3))
  for(k in 1:ncomp){
    print(paste0("####################comp", k))
    beta_inits[[k]] = msCCAproximal_l1$beta_init
    # out = msCCAproximal_l1$direction_update_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R, 
    #                                                beta_init =  msCCAproximal_l1$beta_init, 
    #                                                l1norm_max =  l1norm_max,   l1norm_min = l1norm_min,
    #                                                warm_up =warm_up, trace = T)
    out = msCCAproximal_l1$direction_update_single(X =  msCCAproximal_l1$X, R =  msCCAproximal_l1$R,
                                                   beta_init =  msCCAproximal_l1$beta_init,
                                                   l1norm_max =  l1norm_max,   l1norm_min = l1norm_min,
                                                   warm_up =warm_up, trace = trace)
    errors_track[[k]] = data.frame(matrix(NA, ncol = 3, nrow = ncol(out$beta_augs)))
    errors_track[[k]][,1] = out$beta_norms
    if(step_selection=="penalized"){
      out1 = msCCAproximal_l1$direction_selection(method = "penalized objective", penalty.C = penalty.C, nfolds = nfolds, foldid = NULL, seed = seed, n.core = NULL, multi.core = multi.core)
    }else if(step_selection=="cv"){
      out1 = msCCAproximal_l1$direction_selection(method = "cv", nfolds = nfolds, foldid = foldid, seed = seed, n.core = NULL, multi.core = multi.core)
    }else{
      stop("step_selection must be penalized or cv!")
    }
    errors_track[[k]][,2] = out1$evaluation_obj
    ####
    step_idxs[k] = which.max( out1$evaluation_obj)
    l1bounds[k] = out$bounds[step_idxs[k]]
    if(!is.null(xlist.te)){
      zs = array(NA, dim = c(nrow(xlist.te[[1]]), D, nrow(errors_track[[k]])))
      for(d in 1:D){
        ll = (pss[d]+1):(pss[d+1])
        zs[,d,] = xlist.te[[d]]%*%out$beta_augs[ll,]
      }
      zsum =apply(zs,c(1,3),sum)
      if(k > 1){
        for(j in 1:ncol(zsum)){
          for(k1 in 1:(k-1)){
            zsum[,j] = zsum[,j] - sum(Zsum_orth.te[,k1]*zsum[,j])/sum(Zsum_orth.te[,k1]^2) *Zsum_orth.te[,k1]
          }
        }

      }
      Zsum_orth.te[,k] = zsum[,step_idxs[k]]
      errors_track[[k]][,3] = apply(zsum^2,2,sum)/apply(zs^2,3,sum)
      if(trace){
        print(errors_track[[k]][max(1,step_idxs[k]-25):min(step_idxs[k]+25, nrow(errors_track[[k]])),])
      }
    }
    msCCAproximal_l1$direction_grow(step_idx=step_idxs[k])
    errors_track_selected[k,] = errors_track[[k]][step_idxs[k],]
    if(k < ncomp){
      msCCAproximal_l1$beta_init = msCCAproximal_l1$beta_init_func(msCCAproximal_l1$R)
    }
  }
  
  #calculate Zsum
  return(list(fitted_model = msCCAproximal_l1, beta_inits = beta_inits, beta_norms =l1bounds ,
              step_idxs = step_idxs,  errors_track = errors_track, errors_track_selected = errors_track_selected))
}


riffle_sequential = function(xlist, ncomp, xlist.te = NULL, ss = floor(seq(2, n/4, length.out = 10)), nfolds = 10, foldid = NULL, 
                             n.core = NULL, seed = 2021,  eta = 0.05, maxiter = 5000, multi.core = "lapply"){
  D = length(xlist)
  n = nrow(xlist[[1]])
  ps = sapply(xlist, function(z) dim(z)[2])
  pss = c(0,cumsum(ps))
  rlist = xlist
  ptotal = pss[D+1]
  Lambdahat = matrix(0, nrow = ptotal, ncol = ptotal)
  for(d in 1:D){
    ll = (pss[d]+1):(pss[d+1])
    Lambdahat[ll,ll] = t(xlist[[d]])%*%xlist[[d]]/n
  }
  set.seed(seed)
  if(is.null(foldid)){
    foldid = sample(rep(1:nfolds, each = ceiling(n/nfolds)), n)
  }else{
    nfolds = length(unique(foldid))
  }
  if(is.null(n.core)){
    n.core = min(nfolds, detectCores())
  }
  cv_evaluation = function(fold_id){
    print(paste0("start fold ",fold_id))
    train_id = which(foldid != fold_id)
    ragg_tmp = ragg[train_id,]
    rlist_tmp = list()
    for(i in 1:D){
      rlist_tmp[[i]] = rlist[[i]][train_id,]
    }
    Sigmahat_tmp =t(ragg_tmp)%*%ragg_tmp/n 
    Lambdahat_tmp = matrix(0, nrow = ptotal, ncol = ptotal)
    for(d in 1:D){
      ll = (pss[d]+1):(pss[d+1])
      Lambdahat_tmp[ll,ll] = t(xlist[[d]][train_id,])%*%xlist[[d]][train_id,]/n
    }
    beta_init_tmp = my_init(rlist_tmp, A = 4)
    beta_init_agg =  beta_init_tmp[[1]]
    for(d in 2:D){
      beta_init_agg = c(beta_init_agg, beta_init_tmp[[d]])
    }
    print(paste0("finish initialize fold ",fold_id))
    beta_aggs_tmp = matrix(NA, ncol = length(ss), nrow = ptotal)
    for(l in 1:length(ss)){
      print(l)
      beta_aggs_tmp[,l] = rifle::rifle(A = Sigmahat_tmp, B = Lambdahat_tmp, init = beta_init_agg, k = ss[l],  maxiter =  maxiter)
    }
    test_id = which(foldid == fold_id)
    Zs = array(NA, dim = c(length(test_id), D, length(ss)))
    Zs_residuals = array(NA, dim = c(length(test_id), D, length(ss)))
    for(d in 1:D){
      ll = (pss[d]+1):pss[d+1]
      Zs[,d,] = xlist[[d]][test_id,]%*%beta_aggs_tmp[ll,]
      Zs_residuals[,d,]= rlist[[d]][test_id,]%*%beta_aggs_tmp[ll,]
    }
    Zsums = apply(Zs_residuals,c(1,3),sum)
    B = apply(Zs^2, 3, sum)
    A = apply(Zsums^2, 2, sum)
    rho_components = matrix(0, nrow = length(A), ncol = 2)
    colnames(rho_components) = c("numerator", "denominator")
    rho_components[,1] = A
    rho_components[,2] = B
    print(paste0("end fold id = ", fold_id))
    return(rho_components)
  }
  betas = array(NA, dim = c(ptotal, ncomp))
  errors_track = matrix(NA, nrow = ncomp, ncol = 3)
  Zsum_orth.te = NA
  if(!is.null(xlist.te)){
    Zsum_orth.te = matrix(0, nrow = nrow(xlist.te[[1]]), ncol = ncomp)
  }
  for(k in 1:ncomp){
    print(paste0("####################comp", k))
    ragg = rlist[[1]]
    for(d in 2:D){
      ragg = cbind(ragg, rlist[[d]])
    }
    Sigmahat =t(ragg)%*%ragg/n 
    beta_init = my_init(rlist, A = 4)
    beta_init_agg = beta_init[[1]]
    for(d in 2:D){
      beta_init_agg = c(beta_init_agg, beta_init[[d]])
    }
    beta_aggs = matrix(NA, ncol = length(ss), nrow = ptotal)
    for(l in 1:length(ss)){
      print(paste0("sparsity choice", l))
      beta_aggs[,l] = rifle::rifle(A = Sigmahat, B = Lambdahat, init = beta_init_agg, k = ss[l],  maxiter =  maxiter, eta = eta)
    }
    #################
    #################
    #cv evaluation
    if(multi.core=="mclapply"){
      outputs <-mclapply(1:nfolds,cv_evaluation, mc.cores =n.core)
    }else if(multi.core=="doparallel"){
      outputs <- foreach(i=1:nfolds)%dopar%{
        cv_evaluation(i)
      }
    }else{
      outputs <-lapply(1:nfolds,cv_evaluation)
    }
    evaluation_obj = outputs[[1]]
    for(i in 2:length(outputs)){
      evaluation_obj= evaluation_obj+outputs[[i]]
    }
    rho_hat = evaluation_obj[,1]/evaluation_obj[,2]
    idx_selected = which.max(rho_hat)
    #################
    ##################
    betas[,k] = beta_aggs[,idx_selected]
    errors_track[k,1]=ss[idx_selected]
    errors_track[k,2] = rho_hat[idx_selected]
    #deflation
    if(k==1){
      zsum_out = zsum_func(Xagg = ragg, beta_agg =betas[,k], prev_Zsum_orth = NULL)
      zsum_prev_orth = matrix(zsum_out$zsum_orth, ncol = 1)
      zsum_prev = matrix(zsum_out$zsum, ncol = 1) 
    }else{
      zsum_out = zsum_func(Xagg = ragg, beta_agg =betas[,k], prev_Zsum_orth = zsum_prev_orth)
      zsum_prev = cbind(zsum_prev, zsum_out$zsum)
      zsum_prev_orth = cbind(zsum_prev_orth, zsum_out$zsum_orth)
    }
    deflate_out = deflation_func(Ragg = ragg, R = rlist, zsum_orth = zsum_out$zsum_orth)
    rlist = deflate_out$R
    ragg = deflate_out$Ragg
    #################
    ##################
    #evaluate on test data
    if(!is.null(xlist.te)){
      zs = array(NA, dim = c(nrow(xlist.te[[1]]), D))
      for(d in 1:D){
        ll = (pss[d]+1):(pss[d+1])
        zs[,d] = xlist.te[[d]]%*%betas[ll,k]
      }
      zsum =apply(zs,1,sum)
      if(k > 1){
        for(k1 in 1:(k-1)){
          zsum = zsum - sum(Zsum_orth.te[,k1]*zsum)/sum(Zsum_orth.te[,k1]^2) *Zsum_orth.te[,k1]
        }
      }
      Zsum_orth.te[,k] = zsum
      errors_track[k,3] = sum(zsum^2)/sum(zs^2)
    }
    print(errors_track[k,] )
    #rho_te_eval = diag(t(beta_aggs)%*%A1%*%beta_aggs)/diag(t(beta_aggs)%*%B1%*%beta_aggs)
  }
  return(list(betas = betas, zsum_prev = zsum_prev, zsum_prev_orth= zsum_prev_orth, 
              Zsum_orth.te = Zsum_orth.te, errors_track = errors_track))
  
  
}










