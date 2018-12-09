// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dloglik_likelihood_gradient
double dloglik_likelihood_gradient(int knot, arma::colvec& facility, arma::colvec& delta, arma::mat& z, arma::mat& b_spline, arma::mat& theta, int N);
RcppExport SEXP _SurvTime_dloglik_likelihood_gradient(SEXP knotSEXP, SEXP facilitySEXP, SEXP deltaSEXP, SEXP zSEXP, SEXP b_splineSEXP, SEXP thetaSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type knot(knotSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b_spline(b_splineSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(dloglik_likelihood_gradient(knot, facility, delta, z, b_spline, theta, N));
    return rcpp_result_gen;
END_RCPP
}
// ddloglik_gradient
List ddloglik_gradient(int knot, arma::colvec& facility, arma::mat& z, arma::colvec& delta, arma::mat& b_spline, arma::mat& theta, int number_facility);
RcppExport SEXP _SurvTime_ddloglik_gradient(SEXP knotSEXP, SEXP facilitySEXP, SEXP zSEXP, SEXP deltaSEXP, SEXP b_splineSEXP, SEXP thetaSEXP, SEXP number_facilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type knot(knotSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b_spline(b_splineSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type number_facility(number_facilitySEXP);
    rcpp_result_gen = Rcpp::wrap(ddloglik_gradient(knot, facility, z, delta, b_spline, theta, number_facility));
    return rcpp_result_gen;
END_RCPP
}
// GDboost_gradient
List GDboost_gradient(int knot, double rate, arma::colvec& facility, arma::colvec& delta, arma::mat& z, arma::mat& b_spline, arma::mat theta);
RcppExport SEXP _SurvTime_GDboost_gradient(SEXP knotSEXP, SEXP rateSEXP, SEXP facilitySEXP, SEXP deltaSEXP, SEXP zSEXP, SEXP b_splineSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type knot(knotSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b_spline(b_splineSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(GDboost_gradient(knot, rate, facility, delta, z, b_spline, theta));
    return rcpp_result_gen;
END_RCPP
}
// dloglik_likelihood_stratify
double dloglik_likelihood_stratify(int knot, arma::colvec& facility, arma::colvec& delta, arma::mat& z, arma::mat& b_spline, arma::mat& theta, int N);
RcppExport SEXP _SurvTime_dloglik_likelihood_stratify(SEXP knotSEXP, SEXP facilitySEXP, SEXP deltaSEXP, SEXP zSEXP, SEXP b_splineSEXP, SEXP thetaSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type knot(knotSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b_spline(b_splineSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(dloglik_likelihood_stratify(knot, facility, delta, z, b_spline, theta, N));
    return rcpp_result_gen;
END_RCPP
}
// ddloglik
List ddloglik(int knot, arma::colvec& facility, arma::colvec& delta, arma::mat& z, arma::mat& b_spline, arma::mat& theta, int number_facility);
RcppExport SEXP _SurvTime_ddloglik(SEXP knotSEXP, SEXP facilitySEXP, SEXP deltaSEXP, SEXP zSEXP, SEXP b_splineSEXP, SEXP thetaSEXP, SEXP number_facilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type knot(knotSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b_spline(b_splineSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type number_facility(number_facilitySEXP);
    rcpp_result_gen = Rcpp::wrap(ddloglik(knot, facility, delta, z, b_spline, theta, number_facility));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SurvTime_dloglik_likelihood_gradient", (DL_FUNC) &_SurvTime_dloglik_likelihood_gradient, 7},
    {"_SurvTime_ddloglik_gradient", (DL_FUNC) &_SurvTime_ddloglik_gradient, 7},
    {"_SurvTime_GDboost_gradient", (DL_FUNC) &_SurvTime_GDboost_gradient, 7},
    {"_SurvTime_dloglik_likelihood_stratify", (DL_FUNC) &_SurvTime_dloglik_likelihood_stratify, 7},
    {"_SurvTime_ddloglik", (DL_FUNC) &_SurvTime_ddloglik, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_SurvTime(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
