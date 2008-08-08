
#include "partykit.h"
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

    SEXP prob, index;
    double sum = 0.0;
    int i;

    prob = VECTOR_ELT(split, PROB_SPLIT);
    if (prob != R_NilValue)
        return(prob);
        
    index = index_split(split);
    SET_VECTOR_ELT(split, PROB_SPLIT, prob = allocVector(REALSXP, LENGTH(index)));
    for (i = 0; i < LENGTH(index); i++) {
        REAL(prob)[i] = (double) (INTEGER(index)[i] != NA_INTEGER);
        sum += REAL(prob)[i];
    }
    for (i = 0; i < LENGTH(index); i++)
        REAL(prob)[i] = REAL(prob)[i] / sum;
    return(prob);
}

SEXP info_split(SEXP split) {

    return(VECTOR_ELT(split, INFO_SPLIT));
}

SEXP split_data(SEXP split, SEXP data, SEXP vmatch) {

    if (vmatch == R_NilValue)
        return(VECTOR_ELT(data, varid_split(split) - 1));
    return(VECTOR_ELT(data, INTEGER(vmatch)[varid_split(split) - 1] - 1));
}

double x2d(SEXP x, int obs) {

    double ret = NA_REAL;

    if (isReal(x)) ret = REAL(x)[obs];
    if (isInteger(x)) {
        if (INTEGER(x)[obs] != NA_INTEGER)
            ret = (double) INTEGER(x)[obs];
    }
    /* if (ISNA(ret)) error("can't coerce x to REAL or INTEGER"); */
    return(ret);
}

int kidid_split(SEXP split, SEXP data, SEXP vmatch, int obs) {

    SEXP x, breaks;
    int i, ret = NA_INTEGER;

    x = split_data(split, data, vmatch);
    breaks = breaks_split(split);
    if (breaks == R_NilValue) {
        if (!isInteger(x)) error("x is not an integer");
        ret = INTEGER(x)[obs];
        if (ret != NA_INTEGER) ret = ret - 1;
    } else {
        ret = cut(x2d(x, obs), REAL(breaks), LENGTH(breaks), 
                  right_split(split));
    }

    if (ret != NA_INTEGER) {
        if (index_split(split) != R_NilValue) {
           ret = INTEGER(index_split(split))[ret];
           if (ret != NA_INTEGER) ret = ret - 1;
        }
    }
    return(ret);
}
