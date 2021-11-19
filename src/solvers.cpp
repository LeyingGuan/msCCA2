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

List proximal_gradient_onestep(const List& X,  const List& R, arma::vec& ps, arma::vec&  pss, int ptotal, int n, int D,
                               List& beta, arma::mat& Zs, arma::mat& Zs_residuals, arma::vec& Zsum, 
                               const int norm_type, int l0norm, double l1norm, double l1bound, double l1boundInf,
                               double eta, double eta_ratio, double l1proximal_tol, int l1proximal_maxit, 
                               bool line_search,  int line_maxit, double eta_low, double eps, bool check){
  arma::mat Zs_new = arma::zeros<arma::mat>(n, D);
  arma::mat Zs_residuals_new = arma::zeros<arma::mat>(n, D);
  arma::vec Zsum_new = arma::zeros<arma::vec>(n);
  arma::vec beta_agg = arma::zeros<arma::vec>(ptotal);
  arma::vec beta_agg_new = arma::zeros<arma::vec>(ptotal);
  arma::vec theta_agg = arma::zeros<arma::vec>(ptotal);
  Rcpp::List zeta(D);
  Rcpp::List theta(D);
  Rcpp::List beta_new(D);
  for(int d= 0; d < D; d++){
    arma::vec beta_d = beta[d];
    arma::uvec u_idx = my_range(pss(d),pss(d+1));
    beta_agg(u_idx) = beta_d*1.0;
  }
  double rho = 0.0;
  double A = arma::sum(Zsum%Zsum);
  double B = 0.0;
  for(int d = 0; d < D; d++){
    B += sum(Zs.col(d)%Zs.col(d));
  }
  rho = A/(B+eps);
  for(int d= 0; d < D; d++){
    arma::uvec u_idx = my_range(pss(d),pss(d+1));
    arma::mat Xd = X[d];
    arma::vec beta_d = beta[d];
    arma::vec zeta_d = Xd.t() * Zsum/sqrt(n) - (Xd.t() * Zs.col(d)/sqrt(n)+eps*beta_d) * rho;
    zeta[d] = zeta_d;
  }
  double Bt = 0.0;
  if(norm_type == 0){
    Bt = (l1bound - l1norm);
  }
  double tol =  0.0;
  if(norm_type == 1){
    tol =  0.0;//log(ptotal)/float(n) * Bt * Bt;
  }
  double deterioriation = tol* 2;
  bool update = true;
  double rho_new = 0.0;
  double l1bound_new = 0.0;//change l1bound_old, l1bound_new to lbound_old/lbound_new
  double eta_old = eta*1.2;
  if(eta_old > 0.25){
    eta_old = 0.25;
  }
  double eta_new = eta_old;
  bool check1 = true;
  while(line_search & (deterioriation >= tol) & update & check1){
    check1 = check;
    eta_old = eta_new;
    eta_new = eta_old/1.5;
    if(eta_new < 0.5 * eta_low){
      eta_new = 0.5 * eta_low;
      update = false;
    }
    double Bnew =Bt * (1.0 - eta_old*eta_ratio);
    if(norm_type==1){
      l1bound_new = Bnew+l1norm;
      if(l1bound_new<l1boundInf){
        l1bound_new = l1boundInf;
      }
    }
    for(int d= 0; d < D; d++){
      arma::vec zeta_d = zeta[d];
      arma::vec beta_d = beta[d];
      arma::vec theta_d;
      theta_d =zeta_d * eta_old + beta_d;
      theta_d = zeta_d * eta_old/rho + beta_d;
      arma::uvec u_idx = my_range(pss(d),pss(d+1));
      theta_agg(u_idx) = theta_d;
      beta_agg(u_idx) = beta_d*1.0;
    }
    if(norm_type==1){
      beta_agg_new =l1proximal(theta_agg, l1bound_new,  l1proximal_tol,  l1proximal_maxit);
    }else{
      beta_agg_new =l0proximal(theta_agg, l0norm);
    }
    //Rcout <<   beta_agg_new << "\n";
    Zsum_new =  Zsum_new * 0.0;
    for(int d= 0; d < D; d++){
      arma::uvec u_idx =  my_range(pss(d),pss(d+1));
      arma::vec beta_d = beta_agg_new(u_idx);
      beta_new[d] = beta_d;
      const arma::mat Xd = X[d];
      const arma::mat Rd = R[d];
      Zs_new.col(d) = (Xd * beta_d)/sqrt(n);
      Zs_residuals_new.col(d) = (Rd * beta_d)/sqrt(n);
      Zsum_new +=Zs_residuals_new.col(d);
    }
    A = arma::sum(Zsum_new%Zsum_new);
    B = 0.0;
    for(int d = 0; d < D; d++){
      B += sum(Zs_new.col(d)%Zs_new.col(d));
    }
    rho_new = A/(B+eps);
    deterioriation = (rho-rho_new)/(rho*eta_old*eta_old);
  }
  return List::create(Named("beta") = beta_new,
                      Named("beta_agg") = beta_agg_new,
                      Named("rho") = rho,
                      Named("rho_new") = rho_new,
                      Named("Zs") = Zs_new,
                      Named("Zsum") = Zsum_new,
                      Named("eta") = eta_old,
                      Named("l1bound") = l1bound_new,
                      Named("deterioriation") = deterioriation);
}

// /*
//  We first implement a rank one model.
//  q = 0: l0-constraint
//  q = 1: l1-constraint
//  */
// [[Rcpp::export]]
List msCCA_proximal_rank1(List beta, const List& X,  const List& R, double rho_tol, int rho_maxit,double eta, const int norm_type, int l0norm, 
                    double l1norm, double cl,double eta_ratio, double l1proximal_tol, int l1proximal_maxit,
                    bool line_search,  int line_maxit, double eta_low, double eps, const bool trace = true, const int print_out = 50){
  int D = X.size();
  int n, ptotal = 0;
  arma::vec l1bounds = arma::ones<arma::vec>(rho_maxit);
  l1bounds = l1bounds *  l1norm;
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
  //calculate Zd
  int l0bound = l0norm;
  double tmp =l1norm * sqrt(std::log(ptotal)/n);
  if(tmp > cl){
    tmp = cl;
  }
  l1bounds[0] = l1bounds[0] + std::max(l1norm * tmp, cl *l1norm);
  double l1boundInf = l1norm+l1norm * tmp;
  double err = 10 * rho_tol;
  int it = 0;
  //calculate current Z
  arma::mat Zs = arma::zeros<arma::mat>(n, D);
  arma::mat Zs_residuals = arma::zeros<arma::mat>(n, D);
  arma::vec Zsum = arma::zeros<arma::vec>(n);
  arma::vec beta_agg =  arma::zeros<arma::vec>(ptotal);
  for(int d= 0; d < D; d++){
    arma::vec beta_d = beta[d];
    //arma::uvec u_idx =  feature_idx[d];
    arma::uvec u_idx =  my_range(pss(d),pss(d+1));
    beta_agg(u_idx) = beta_d*1.0;
  }
  for(int d = 0; d < D; d++){
    arma::vec beta_d = beta[d];
    const arma::mat Xd = X[d];
    const arma::mat Rd = R[d];
    Zs.col(d) = (Xd * beta_d)/sqrt(n);
    Zs_residuals.col(d) = (Rd * beta_d)/sqrt(n);
    Zsum +=Zs_residuals.col(d);
  }
  double rho = 0.0;
  double rho_new = 0.0;
  double rho_change = 0.0;
  double l1bound = 0.0;
  bool check = false;
  while(err > rho_tol & it < (rho_maxit-1)){
    l1bound = l1bounds[it];
    if(it>0){
      check = true;
    }
    List out = proximal_gradient_onestep(X, R, ps, pss, ptotal, n, D, beta, Zs, Zs_residuals, Zsum,
                                         norm_type, l0norm, l1norm, l1bound, l1boundInf,
                                         eta, eta_ratio, l1proximal_tol, l1proximal_maxit,
                                         line_search, line_maxit, eta_low,  eps, check);
    eta = out["eta"];
    arma::mat Zs_new =  out["Zs"];
    arma::vec Zsum_new = out["Zsum"];
    rho = out["rho"];
    rho_new = out["rho_new"];
    arma::vec beta_agg_new = out["beta_agg"];
    double l1bound_new = out["l1bound"];
    Zs = Zs_new * 1.0;
    Zsum = Zsum_new * 1.0;
    for(int d = 0; d < D; d++){
      arma::uvec u_idx =  my_range(pss(d),pss(d+1));
      arma::vec beta_d = beta_agg_new(u_idx);
      arma::vec beta_d_copy = beta_d *1.0;
      beta[d] = beta_d_copy;
    }
    arma::vec change = beta_agg_new - beta_agg;
    err = arma::sum(change%change)/(eta*eta);
    beta_agg = beta_agg_new * 1.0;
    rho_change = (rho - rho_new)/rho;
    if(trace){
      if(it%print_out == 0){
        Rcout << "rho_iteration: " << it << ", eta: " <<  eta << ", new obj: " <<rho_new << ", relative beta change (beta/eta):" << err << "\n";
      }
    }
    it += 1;
    l1bounds[it] = l1bound_new;
  }
  l1bound = l1bounds[it];
  double a = arma::sum(Zsum%Zsum);
  double b = 0.0;
  for(int d = 0; d < D; d++){
    b+= arma::sum(Zs.col(d)%Zs.col(d));
  }
  rho_new = a/b;
  List tmp_return =  List::create(Named("beta") = beta,
                                  Named("beta_aug") = beta_agg,
                                  Named("Zs") = Zs,
                                  Named("Zsum") = Zsum,
                                  Named("rho") = rho_new,
                                  Named("eta") = eta,
                                  Named("l1bound") = l1bound);
  return tmp_return;
}
