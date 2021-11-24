#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _msCCA2_msCCA_proximal_rank1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _msCCA2_my_range(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_msCCA2_msCCA_proximal_rank1", (DL_FUNC) &_msCCA2_msCCA_proximal_rank1, 20},
    {"_msCCA2_my_range",             (DL_FUNC) &_msCCA2_my_range,              2},
    {NULL, NULL, 0}
};

void R_init_msCCA2(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}