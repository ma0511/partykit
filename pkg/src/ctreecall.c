
#include "ctree.h"

SEXP R_LinstatExpCov (const SEXP x, const SEXP y, const SEXP weights) {

    SEXP ans;
    int *thisweights;

    thisweights = Calloc(LENGTH(x), int);
    PROTECT(ans = allocVector(VECSXP, 4));
    C_LinstatExpCov(x, y, weights, thisweights, ans);
    UNPROTECT(1);
    Free(thisweights);
    return(ans);
}

SEXP R_split (const SEXP x, const SEXP y, const SEXP weights, const SEXP minbucket) {

    SEXP ans, cx;
    int i, n, q, sw = 0;
    int *thisweights;
    double *expinf, *covinf;

    n = NROW(x);
    q = NCOL(y);     

    expinf = Calloc(q, double);
    covinf = Calloc(q * q, double);
    thisweights = Calloc(n, int);
    
    if (isReal(x)) {
        PROTECT(ans = allocVector(REALSXP, 1));
        NA_weights_double(INTEGER(weights), REAL(x), n, thisweights, &sw);
        C_split_numeric(x, y, thisweights, INTEGER(minbucket)[0], expinf, covinf, REAL(ans));
    }
    if (isFactor(x) && !isOrdered(x)) {
        PROTECT(ans = allocVector(INTSXP, C_nlevels(x)));
        NA_weights_factor(INTEGER(weights), INTEGER(x), n, thisweights, &sw);
        C_split_factor(x, y, thisweights, INTEGER(minbucket)[0], expinf, covinf, INTEGER(ans));
    }
    if (isInteger(x) || isOrdered(x)) {
        PROTECT(ans = allocVector(REALSXP, 1));
        PROTECT(cx = coerceVector(x, REALSXP));
        NA_weights_double(INTEGER(weights), REAL(cx), n, thisweights, &sw);
        C_split_numeric(cx, y, thisweights, INTEGER(minbucket)[0], expinf, covinf, REAL(ans));
        UNPROTECT(1);
    }

    Free(thisweights);
    Free(expinf);
    Free(covinf);

    UNPROTECT(1);
    return(ans);
}
