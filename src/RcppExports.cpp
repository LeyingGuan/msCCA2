// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// my_range
arma::uvec my_range(int n1, int n2);
RcppExport SEXP _msCCA_my_range(SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(my_range(n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// msCCA_proximal_rank1
List msCCA_proximal_rank1(List beta, const List& X, const List& R, double rho_tol, int rho_maxit, double eta, const int norm_type, int l0norm, double l1norm, double cl, double eta_ratio, double l1proximal_tol, int l1proximal_maxit, bool line_search, int line_maxit, double eta_low, double eps, const bool trace, const int print_out);
RcppExport SEXP _msCCA_msCCA_proximal_rank1(SEXP betaSEXP, SEXP XSEXP, SEXP RSEXP, SEXP rho_tolSEXP, SEXP rho_maxitSEXP, SEXP etaSEXP, SEXP norm_typeSEXP, SEXP l0normSEXP, SEXP l1normSEXP, SEXP clSEXP, SEXP eta_ratioSEXP, SEXP l1proximal_tolSEXP, SEXP l1proximal_maxitSEXP, SEXP line_searchSEXP, SEXP line_maxitSEXP, SEXP eta_lowSEXP, SEXP epsSEXP, SEXP traceSEXP, SEXP print_outSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const List& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const List& >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type rho_tol(rho_tolSEXP);
    Rcpp::traits::input_parameter< int >::type rho_maxit(rho_maxitSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const int >::type norm_type(norm_typeSEXP);
    Rcpp::traits::input_parameter< int >::type l0norm(l0normSEXP);
    Rcpp::traits::input_parameter< double >::type l1norm(l1normSEXP);
    Rcpp::traits::input_parameter< double >::type cl(clSEXP);
    Rcpp::traits::input_parameter< double >::type eta_ratio(eta_ratioSEXP);
    Rcpp::traits::input_parameter< double >::type l1proximal_tol(l1proximal_tolSEXP);
    Rcpp::traits::input_parameter< int >::type l1proximal_maxit(l1proximal_maxitSEXP);
    Rcpp::traits::input_parameter< bool >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< int >::type line_maxit(line_maxitSEXP);
    Rcpp::traits::input_parameter< double >::type eta_low(eta_lowSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const int >::type print_out(print_outSEXP);
    rcpp_result_gen = Rcpp::wrap(msCCA_proximal_rank1(beta, X, R, rho_tol, rho_maxit, eta, norm_type, l0norm, l1norm, cl, eta_ratio, l1proximal_tol, l1proximal_maxit, line_search, line_maxit, eta_low, eps, trace, print_out));
    return rcpp_result_gen;
END_RCPP
}