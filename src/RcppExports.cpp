// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sdRcpp
NumericVector sdRcpp(NumericMatrix x);
RcppExport SEXP _ddsPLS_sdRcpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sdRcpp(x));
    return rcpp_result_gen;
END_RCPP
}
// get_sd_matrixRcpp
NumericMatrix get_sd_matrixRcpp(NumericMatrix x);
RcppExport SEXP _ddsPLS_get_sd_matrixRcpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sd_matrixRcpp(x));
    return rcpp_result_gen;
END_RCPP
}
// scaleRcpp
NumericMatrix scaleRcpp(NumericMatrix x);
RcppExport SEXP _ddsPLS_scaleRcpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(scaleRcpp(x));
    return rcpp_result_gen;
END_RCPP
}
// mmultC
NumericMatrix mmultC(NumericMatrix m1, NumericMatrix m2);
RcppExport SEXP _ddsPLS_mmultC(SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(mmultC(m1, m2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ddsPLS_sdRcpp", (DL_FUNC) &_ddsPLS_sdRcpp, 1},
    {"_ddsPLS_get_sd_matrixRcpp", (DL_FUNC) &_ddsPLS_get_sd_matrixRcpp, 1},
    {"_ddsPLS_scaleRcpp", (DL_FUNC) &_ddsPLS_scaleRcpp, 1},
    {"_ddsPLS_mmultC", (DL_FUNC) &_ddsPLS_mmultC, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ddsPLS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
