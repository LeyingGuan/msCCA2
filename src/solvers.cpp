#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//return ranks from large to small
// [[Rcpp::export]]
arma::uvec my_range(int n1, int n2){
  arma::uvec u_idx = arma::zeros<arma::uvec>(n2-n1);
  for(int i = n1; i < n2; i++){
    u_idx[i-n1] = i;
  }
  return u_idx;
}

List get_abs_ranks(arma::vec theta_abs){
  arma::uvec sorted_index = arma::sort_index(-theta_abs);
  arma::vec  theta_abs_sorted =theta_abs(sorted_index);
  arma::uvec reverse_ranks = arma::sort_index(sorted_index);
  return List::create(Named("theta_abs_sorted") = theta_abs_sorted,
                      Named("sorted_index") = sorted_index);
}


arma::vec l1proximal(arma::vec theta, double l1bound, double tol, int max_iter){
  //get the rank of |theta|
  arma::vec theta_abs = abs(theta);
  arma::vec theta_signs = sign(theta);
  List sort_info = get_abs_ranks(theta_abs);
  arma::vec theta_abs_sorted = sort_info["theta_abs_sorted"];
  arma::uvec sorted_index = sort_info["sorted_index"];
  int p = theta.n_elem;
  int a1 =  floor(l1bound * l1bound)-1;
  int s = std::min(a1, p-1);
  double lower_l1 = 0.0;
  double upper_l1 = std::min(theta_abs_sorted[s],theta_abs_sorted[0]-1e-8);
  double middle = (upper_l1 + lower_l1)/2.0;
  //check if unpenalized version works
  arma::vec theta_abs_thrs = theta_abs;
  theta_abs_thrs =theta_abs_thrs/sqrt(arma::sum(theta_abs_thrs%theta_abs_thrs));
  double l1realized = sum(theta_abs_thrs);
  if(l1realized > l1bound){
    theta_abs_thrs = (theta_abs - middle);
    arma::uvec tmp_idx = arma::find(theta_abs_thrs < 0.0);
    theta_abs_thrs(tmp_idx) = theta_abs_thrs(tmp_idx) * 0.0;
    theta_abs_thrs =theta_abs_thrs/sqrt(arma::sum(theta_abs_thrs%theta_abs_thrs));
    l1realized = sum(theta_abs_thrs);
    double err = l1realized - l1bound;
    int it = 0;
    while(((err > tol) | (err < -tol)) & (it < max_iter)){
      it += 1;
      if(err > tol){
        //realized l1 norm too large
        lower_l1 = middle;
        middle = (upper_l1+ lower_l1)/2.0;
      }else if(err < -tol){
        //realized l1 norm too small
        upper_l1 = middle;
        middle = (upper_l1+ lower_l1)/2.0;
      }
      theta_abs_thrs = (theta_abs - middle);
      tmp_idx = arma::find(theta_abs_thrs < 0.0);
      theta_abs_thrs(tmp_idx) = theta_abs_thrs(tmp_idx) * 0.0;
      theta_abs_thrs =theta_abs_thrs/sqrt(arma::sum(theta_abs_thrs%theta_abs_thrs));
      l1realized = sum(theta_abs_thrs);
      err = l1realized - l1bound;
    }
    if(abs(err) > tol){
      Rcout << "binary search does not produce desired precision, and more iterations might be needed." << "\n";
    }
  }
  arma::vec theta_thrs = theta_abs_thrs % theta_signs;
  return theta_thrs;
}

arma::vec l0proximal(arma::vec theta, int l0bound){
  //keep the l0bound entries with the largest absolute magnitute and zeros out others
  arma::vec theta_abs = abs(theta);
  arma::vec theta_signs = sign(theta);
  List sort_info = get_abs_ranks(theta_abs);
  arma::vec theta_abs_sorted = sort_info["theta_abs_sorted"];
  arma::uvec sorted_index = sort_info["sorted_index"];
  int p = theta.n_elem;
  double thr = theta_abs_sorted[l0bound];
  for(int i =0; i < p; i++){
    if(theta_abs[i]<thr){
      theta[i] = 0.0;
    }
  }
  theta = theta/sqrt(arma::sum(theta%theta));
  return theta;
}

List evaluation(arma::vec& beta_agg, const List& X,  const List& R,  int ptotal, int n, int D,
                arma::vec& ps, arma::vec&  pss, double eps){
  arma::mat Zs = arma::zeros<arma::mat>(n, D);
  arma::mat Zs_residuals = arma::zeros<arma::mat>(n, D);
  arma::vec Zsum = arma::zeros<arma::vec>(n);
  arma::vec zeta_agg = arma::zeros<arma::vec>(ptotal);
  Rcpp::List beta(D);
  for(int d= 0; d < D; d++){
    arma::uvec u_idx =  my_range(pss(d),pss(d+1));
    arma::vec beta_d = beta_agg(u_idx);
    beta[d] = beta_d;
    const arma::mat Xd = X[d];
    const arma::mat Rd = R[d];
    Zs.col(d) = (Xd * beta_d)/sqrt(n);
    Zs_residuals.col(d) = (Rd * beta_d)/sqrt(n);
    Zsum +=Zs_residuals.col(d);
  }
  double A = arma::sum(Zsum%Zsum);
  double B = 0.0;
  for(int d = 0; d < D; d++){
    B += sum(Zs.col(d)%Zs.col(d));
  }
  double rho = A/(B+eps);
  for(int d = 0; d < D; d++){
    arma::uvec u_idx =  my_range(pss(d),pss(d+1));
    arma::vec beta_d = beta_agg(u_idx);
    const arma::mat Xd = X[d];
    const arma::mat Rd = R[d];
    zeta_agg(u_idx) = Rd.t() * Zsum/sqrt(n) - (Xd.t()*Zs.col(d)/sqrt(n)+eps * beta_d)*rho;
  }
  return List::create(Named("beta") = beta,
                      Named("rho_denominator") = B,
                      Named("rho_numerator") = A,
                      Named("zeta_agg") = zeta_agg);
}

List proximal_gradient_onestep(const List& X,  const List& R, arma::vec& ps, arma::vec&  pss, int ptotal, int n, int D,
                               arma::vec& beta_agg, arma::vec& zeta_agg, double rho_prev,
                               const int norm_type, int l0norm, double l1norm_min, double l1bound, 
                               double eta, double eta_ratio, double l1proximal_tol, int l1proximal_maxit, int line_maxit, 
                               double eta_low, double eps){
  arma::vec beta_agg_new = arma::zeros<arma::vec>(ptotal);
  arma::vec theta_agg = arma::zeros<arma::vec>(ptotal);
  double A, B;
  double tol =  -1e-8;
  // if(norm_type == 1){
  //   tol =  0.0;//log(ptotal)/float(n) * Bt * Bt;
  // }
  double rho_new = 0.0;
  double bound_new = 0.0;//change l1bound_old, l1bound_new to lbound_old/lbound_new
  //proximal update and decide if changing step size:
  // if increase step size leads to more objective improvement, increase 20%;
  // if the objective deterioate, deacrease 50%
  //decide the step size using old bound
  theta_agg = beta_agg + zeta_agg * eta/rho_prev;
  if(norm_type == 1){
    bound_new = (l1bound - l1norm_min) * (1.0 - eta* eta_ratio)+l1norm_min;
    beta_agg_new =l1proximal(theta_agg, bound_new,  l1proximal_tol,  l1proximal_maxit);
  }else{
    bound_new = l0norm;
    beta_agg_new = l0proximal(theta_agg, bound_new);
  }
  List eval_out = evaluation(beta_agg_new, X,  R,  ptotal, n, D, ps, pss, eps);
  A = eval_out["rho_numerator"];
  B = eval_out["rho_denominator"];
  rho_new = A/(B+eps);
  arma::vec zeta_agg_new = eval_out["zeta_agg"];
  //check if the current step size eta is too large:
  double improvement = (rho_new - rho_prev);
  // bool shrink = true;
  // if(improvement >= tol){
  //   shrink = false;
  // }
  // if(shrink){
  //   //intrigue line search
  //   bool update = true;
  //   if(eta < eta_low){
  //     update = false;
  //   }
  //   while(update){
  //     eta = eta/2.0;
  //     theta_agg = beta_agg + zeta_agg * eta/rho_prev;
  //     if(norm_type == 1){
  //       bound_new = l1bound;
  //       beta_agg_new =l1proximal(theta_agg, bound_new,  l1proximal_tol,  l1proximal_maxit);
  //     }else{
  //       bound_new = l0norm;
  //       beta_agg_new = l0proximal(theta_agg, bound_new);
  //     }
  //     eval_out = evaluation(beta_agg_new, X,  R,  ptotal, n, D, ps, pss, eps);
  //     A = eval_out["rho_numerator"];
  //     B = eval_out["rho_denominator"];
  //     rho_new = A/(B+eps);
  //     improvement = (rho_new - rho_prev);
  //     if(eta < eta_low){
  //       update = false;
  //     }else if(improvement > tol){
  //       update = false;
  //     }
  //   }
  // }
  // arma::vec beta_agg_new1 = arma::zeros<arma::vec> (ptotal);
  // double A1, B1, bound_new1;
  // arma::vec theta_agg1 = beta_agg + zeta_agg * eta*1.2/rho_prev;
  // if(norm_type==1){
  //   bound_new1 = l1bound;
  //   beta_agg_new1 =l1proximal(theta_agg1, bound_new1,  l1proximal_tol,  l1proximal_maxit);
  // }else{
  //   bound_new1 = l0norm;
  //   beta_agg_new1 = l0proximal(theta_agg1, bound_new1);
  // }
  // List eval_out1 = evaluation(beta_agg_new1, X,  R,  ptotal, n, D, ps, pss, eps);
  // A1 = eval_out["rho_numerator"];
  // B1 = eval_out["rho_denominator"];
  // double rho_new1 = A1/(B1+eps);
  // if(rho_new1 > rho_new){
  //   eta = eta* 1.2;
  // }
  // step 1: update with current step size
  // theta_agg = beta_agg + zeta_agg * eta/rho_prev;
  // if(norm_type == 1){
  //   bound_new = (l1bound - l1norm_min) * (1.0 - eta* eta_ratio)+l1norm_min;
  //   beta_agg_new =l1proximal(theta_agg, bound_new,  l1proximal_tol,  l1proximal_maxit);
  // }else{
  //   bound_new = l0norm;
  //   beta_agg_new = l0proximal(theta_agg, bound_new);
  // }
  // eval_out = evaluation(beta_agg_new, X,  R,  ptotal, n, D, ps, pss, eps);
  // A = eval_out["rho_numerator"];
  // B = eval_out["rho_denominator"];
  // arma::vec zeta_agg_new = eval_out["zeta_agg"];
  // rho_new = A/(B+eps);
  // improvement = (rho_new - rho_prev);
  return List::create(Named("beta_agg") = beta_agg_new,
                      Named("rho") = rho_new,
                      Named("rho_denominator") = B,
                      Named("rho_numerator") = A,
                      Named("zeta_agg") = zeta_agg_new,
                      Named("eta") = eta,
                      Named("bound") = bound_new,
                      Named("improvement") = improvement);
}

// /*
//  We first implement a rank one model.
//  q = 0: l0-constraint
//  q = 1: l1-constraint
//  new functions:
//  track the achieved empirical rho on training
//  track the l1 norm
//  track the achieved numerator and denominator on test
//  */
// [[Rcpp::export]]
List msCCA_proximal_rank1(List beta, const List& X,  const List& R, double rho_tol, int rho_maxit,double eta, 
                          const int norm_type, int l0norm, double  l1norm_max, double l1norm_min, 
                           double eta_ratio, double l1proximal_tol, int l1proximal_maxit,
                     int line_maxit, double eta_low, double eps, int warm_up,
                    const bool trace = true, const int print_out = 50, bool early_stop = true){
  int D = X.size();
  int n, ptotal = 0;
  arma::vec ps = arma::zeros<arma::vec>(D);
  arma::vec pss = arma::zeros<arma::vec>(D+1);
  for(int d = 0; d < D; d++){
    const arma::mat Xd = X[d];
    arma::vec beta_d = beta[d];
    ps(d) = Xd.n_cols;
    n = Xd.n_rows;
    ptotal += ps(d);
    pss(d+1) = pss(d)+ps(d);
  }
  arma::vec beta_agg =  arma::zeros<arma::vec>(ptotal);
  for(int d= 0; d < D; d++){
    arma::vec beta_d = beta[d];
    //arma::uvec u_idx =  feature_idx[d];
    arma::uvec u_idx =  my_range(pss(d),pss(d+1));
    beta_agg(u_idx) = beta_d*1.0;
  }
  arma::vec beta_norms = arma::ones<arma::vec>(rho_maxit);
  arma::vec train_rhos = arma::ones<arma::vec>(rho_maxit);
  arma::vec rho_denominators = arma::ones<arma::vec>(rho_maxit);
  arma::vec rho_numerators = arma::ones<arma::vec>(rho_maxit);
  arma::vec bounds = arma::ones<arma::vec>(rho_maxit);
  arma::vec etas = arma::ones<arma::vec>(rho_maxit);
  arma::mat betas_all = arma::zeros<arma::mat>(ptotal, rho_maxit);
  if(norm_type == 1){
    bounds = bounds * l1norm_max;
  }else{
    bounds = bounds * l0norm;
  }
  int warm_up_it = 0;
  int it = 0;
  beta_norms(it) = sum(abs(beta_agg));
  betas_all.col(it) =beta_agg*1.0; 
  double rho_denominator, rho_numerator;
  List eval_out = evaluation(beta_agg, X,  R,  ptotal, n, D, ps,  pss,  eps);
  rho_denominator = eval_out["rho_denominator"];
  rho_numerator = eval_out["rho_numerator"];
  arma::vec zeta_agg = eval_out["zeta_agg"];
  rho_denominators(it) = rho_denominator;
  rho_numerators(it) = rho_numerator;
  train_rhos(it) =  rho_numerator/(rho_denominator+eps);
  double rho = train_rhos(it);
  double bound = bounds(it);
  while(warm_up_it < warm_up){
    List out = proximal_gradient_onestep(X, R, ps, pss, ptotal, n, D, beta_agg,  zeta_agg, rho, 
                                         norm_type, l0norm, l1norm_min, bound,
                                         eta, 0, l1proximal_tol, l1proximal_maxit,
                                         line_maxit, eta_low,  eps);
    arma::vec beta_agg_new = out["beta_agg"] ;
    arma::vec zeta_agg_tmp = out["zeta_agg"] ;
    zeta_agg= zeta_agg_tmp*1.0;
    beta_agg = beta_agg_new * 1.0;
    rho_denominator = out["rho_denominator"];
    rho_numerator = out["rho_numerator"];
    rho = rho_numerator/(rho_denominator+eps);
    warm_up_it+=1;
  }
  rho_denominators(it) = rho_denominator;
  rho_numerators(it) = rho_numerator;
  train_rhos(it) =  rho_numerator/(rho_denominator+eps);
  beta_norms(it) = sum(abs(beta_agg));
  betas_all.col(it) =beta_agg*1.0; 
  double relative_change = 10 * rho_tol;
  double bound_new;
  etas(it) = eta*1.0;
  while(relative_change > rho_tol & it < (rho_maxit-1)){
    bound = bounds[it] * 1.0;
    rho = train_rhos(it) * 1.0;
    beta_agg = betas_all.col(it) * 1.0;
    List out = proximal_gradient_onestep(X, R, ps, pss, ptotal, n, D, beta_agg,  zeta_agg, rho, 
                                         norm_type, l0norm, l1norm_min, bound,
                                         eta, eta_ratio, l1proximal_tol, l1proximal_maxit,
                                         line_maxit, eta_low,  eps);
    arma::vec beta_agg_new = out["beta_agg"] ;
    arma::vec zeta_agg_tmp = out["zeta_agg"] ;
    zeta_agg= zeta_agg_tmp*1.0;
    bound_new = out["bound"] ;
    rho_denominator = out["rho_denominator"];
    rho_numerator = out["rho_numerator"];
    eta = out["eta"];
    it += 1;
    betas_all.col(it) = beta_agg_new * 1.0;
    beta_norms(it) = arma::sum(arma::abs(beta_agg_new));
    bounds(it) = bound_new;
    rho_denominators(it) = rho_denominator;
    rho_numerators(it) = rho_numerator;
    train_rhos(it) =  rho_numerator/(rho_denominator+eps);
    etas(it) = eta*1.0;
    double beta_change = sqrt(arma::sum(arma::square(betas_all.col(it) - betas_all.col(it-1))));
    relative_change = beta_change/eta;
    if(trace){
      if(it%print_out == 0){
        Rcout << "rho_iteration: " << it << ", eta: " <<  eta << ", new obj: " <<train_rhos(it) << ", relative beta change (beta/eta):" << relative_change<< "\n";
      }
    }
    if(!early_stop){
      relative_change = rho_tol * 10;
    }
  }
  List tmp_return =  List::create(Named("beta_augs") = betas_all,
                                  Named("bounds") = bounds,
                                  Named("beta_norms") = beta_norms,
                                  Named("rhos") = train_rhos,
                                  Named("rho_denominators")=rho_denominators,
                                  Named("rho_numerators")=rho_numerators,
                                  Named("etas") = etas,
                                  Named("end_iter") = it
                                  );
  return tmp_return;
}


