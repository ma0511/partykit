
#include "split.h"

void init_split(SEXP varid, SEXP breaks, SEXP index, SEXP right,
             SEXP prob, SEXP info, SEXP split) {
             
    if (LENGTH(split) != LENGTH_SPLIT)
        error("split is not a list with %d elements", LENGTH_SPLIT);
        
    if (isInteger(varid))
        SET_VECTOR_ELT(split, VARID_SPLIT, varid);

    /* FIXME: sort first? */
    SET_VECTOR_ELT(split, BREAKS_SPLIT, coerceVector(breaks, REALSXP));

    SET_VECTOR_ELT(split, INDEX_SPLIT, index);
    SET_VECTOR_ELT(split, RIGHT_SPLIT, right);
    SET_VECTOR_ELT(split, PROB_SPLIT, prob);
    SET_VECTOR_ELT(split, INFO_SPLIT, info);
}

SEXP R_split(SEXP varid, SEXP breaks, SEXP index, SEXP right,
             SEXP prob, SEXP info) {
             
    SEXP split;
    
    PROTECT(split = allocVector(VECSXP, LENGTH_SPLIT));
    init_split(varid, breaks, index, right, prob, info, split);
    UNPROTECT(1);
    return(split);
}

int varid_split(SEXP split) {

    return(INTEGER(VECTOR_ELT(split, VARID_SPLIT))[0]);
}

SEXP breaks_split(SEXP split) {

    return(VECTOR_ELT(split, BREAKS_SPLIT));
}

SEXP index_split(SEXP split) {

    return(VECTOR_ELT(split, INDEX_SPLIT));
}

int right_split(SEXP split) {

    return(LOGICAL(VECTOR_ELT(split, RIGHT_SPLIT))[0]);
}

SEXP prob_split(SEXP split) {

    return(VECTOR_ELT(split, PROB_SPLIT));
}

SEXP info_split(SEXP split) {

    return(VECTOR_ELT(split, INFO_SPLIT));
}

SEXP split_data(SEXP split, SEXP data, SEXP vmatch) {

    if (vmatch == R_NilValue)
        return(VECTOR_ELT(data, varid_split(split) - 1));
    return(VECTOR_ELT(data, INTEGER(vmatch)[varid_split(split) - 1] - 1));
}

int cut(double x, double *breaks, int n, int right) {

    int ret, i;

    ret = NA_INTEGER;
    if (x > breaks[n - 1]) {
        ret = n;
    } else {
        for (i = 0; i < n; i++) {
            if (x <= breaks[i]) {
                ret = i;
                break;
            }
        }
        if (!right)
            if (x == breaks[ret]) ret++;
    }
    return(ret);
}

double x2d(SEXP x, int obs) {

    double ret = NA_REAL;
    if (isReal(x)) ret = REAL(x)[obs];
    if (isInteger(x)) ret = (double) INTEGER(x)[obs];
    /* if (ISNA(ret)) error("can't coerce x to REAL or INTEGER"); */
    return(ret);
}

int kidid_split(SEXP split, SEXP data, SEXP vmatch, int obs) {

    SEXP x, breaks;
    int i, ret = NA_INTEGER;
    double dx;

    x = split_data(split, data, vmatch);
    breaks = breaks_split(split);
    if (breaks == R_NilValue) {
        if (!isInteger(x)) error("x is not an integer");
        ret = INTEGER(x)[obs] - 1;
    } else {
        ret = cut(x2d(x, obs), REAL(breaks), LENGTH(breaks), 
                  right_split(split));
    }

    if (ret != NA_INTEGER) {
        if (index_split(split) != R_NilValue)
           ret = INTEGER(index_split(split))[ret] - 1;
    }

    return(ret);
}

SEXP R_kidids_split(SEXP split, SEXP data, SEXP vmatch, SEXP obs) {

    SEXP ans;
    int i, tmp;   
        
    PROTECT(ans = allocVector(INTSXP, LENGTH(obs)));
    for (i = 0; i < LENGTH(ans); i++) {
        tmp = kidid_split(split, data, vmatch, INTEGER(obs)[i] - 1);
        if (tmp != NA_INTEGER) {
            INTEGER(ans)[i] = tmp + 1;
        } else {
            INTEGER(ans)[i] = NA_INTEGER;
        }
        
    }
    UNPROTECT(1);
    return(ans);
}
