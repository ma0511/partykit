
/**
    .Call interfaces
    *\file utils.c
    *\author $Author$
    *\date $Date$
*/

#include "party.h"
#include "partysplit.h"
#include "partynode.h"

/**
    determine the child node observations obs has to go into based
    on one split (either primary or surrogate)
    *\param split a split object
    *\param data a list
    *\param vmatch an integer for permuting variables
    *\param obs integer vector of observation numbers
*/
                
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

/**
    determine the terminal node id for observations obs.
    *\param node a node object
    *\param data a list
    *\param vmatch an integer for permuting variables
    *\param obs integer vector of observation numbers
*/

SEXP R_fitted_node(SEXP node, SEXP data, SEXP vmatch, SEXP obs) {

    SEXP ans;
    int i;

    /* we might want to do random splitting */
    GetRNGstate();
         
    PROTECT(ans = allocVector(INTSXP, LENGTH(obs)));
    for (i = 0; i < LENGTH(ans); i++)
        INTEGER(ans)[i] = fitted_node(node, data, vmatch, INTEGER(obs)[i] - 1);

    PutRNGstate();

    UNPROTECT(1);
    return(ans);
}
