/***************************************************************************\
 *
 * Authors: Guillaume Guénard
 * Université de Montréal - 2023
 *
 * Simulation of DNA sequences evolution as a Markov process at the scale of
 * single nucleotides, with the purpose of generating simulated data for deep
 * learning model developments.
 *
 * Dynamic symbols exports
 *
 \***************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dstIdx_Cpp
IntegerVector dstIdx_Cpp(const int, const IntegerVector&,
                         const IntegerVector&);
RcppExport SEXP PhyloSim_dstIdx(SEXP nSEXP, SEXP fromSEXP, SEXP toSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< const int >::type n(nSEXP);
  Rcpp::traits::input_parameter< const IntegerVector& >::type from(fromSEXP);
  Rcpp::traits::input_parameter< const IntegerVector& >::type to(toSEXP);
  rcpp_result_gen = Rcpp::wrap(dstIdx_Cpp(n, from, to));
  return rcpp_result_gen;
  END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
  {"PhyloSim_dstIdx", (DL_FUNC) &PhyloSim_dstIdx, 3},
  {NULL, NULL, 0}
};

RcppExport void R_init_PhyloSim(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
