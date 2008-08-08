
#include "partykit.h"
#include "split.h"
#include "node.h"

SEXP R_split(SEXP varid, SEXP breaks, SEXP index, SEXP right,
             SEXP prob, SEXP info) {
             
    SEXP split;
    
    PROTECT(split = allocVector(VECSXP, LENGTH_SPLIT));
    init_split(varid, breaks, index, right, prob, info, split);
    UNPROTECT(1);
    return(split);
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

SEXP R_fitted_node(SEXP node, SEXP data, SEXP vmatch, SEXP obs) {

    SEXP ans;
    int i;

    PROTECT(ans = allocVector(INTSXP, LENGTH(obs)));
    for (i = 0; i < LENGTH(ans); i++)
        INTEGER(ans)[i] = fitted_node(node, data, vmatch, INTEGER(obs)[i] - 1);
    UNPROTECT(1);
    return(ans);
}
