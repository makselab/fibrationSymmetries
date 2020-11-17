// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getBalancedColoring
std::vector<int> getBalancedColoring(std::vector<int> nodeColors, std::vector<bool> colorFixed, IntegerMatrix edges, bool directed, bool weighted, int numberOfWeights);
RcppExport SEXP _fibrationSymmetries_getBalancedColoring(SEXP nodeColorsSEXP, SEXP colorFixedSEXP, SEXP edgesSEXP, SEXP directedSEXP, SEXP weightedSEXP, SEXP numberOfWeightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type nodeColors(nodeColorsSEXP);
    Rcpp::traits::input_parameter< std::vector<bool> >::type colorFixed(colorFixedSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type edges(edgesSEXP);
    Rcpp::traits::input_parameter< bool >::type directed(directedSEXP);
    Rcpp::traits::input_parameter< bool >::type weighted(weightedSEXP);
    Rcpp::traits::input_parameter< int >::type numberOfWeights(numberOfWeightsSEXP);
    rcpp_result_gen = Rcpp::wrap(getBalancedColoring(nodeColors, colorFixed, edges, directed, weighted, numberOfWeights));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fibrationSymmetries_getBalancedColoring", (DL_FUNC) &_fibrationSymmetries_getBalancedColoring, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_fibrationSymmetries(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}