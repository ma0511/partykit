
#include "partykit.h"
#include "node.h"
#include "split.h"

/*
node <- function(id, split = NULL, kids = NULL, surrogates = NULL, info = NULL) {
    node <- list(id = id, split = split, kids = kids, surrogates = surrogates, 	info = info)
}
*/

int id_node(SEXP node) {
    return(INTEGER(VECTOR_ELT(node, ID_NODE))[0]);
}

SEXP split_node(SEXP node) {
    return(VECTOR_ELT(node, SPLIT_NODE));
}

SEXP kids_node(SEXP node) {
    return(VECTOR_ELT(node, KIDS_NODE));
}

SEXP surrogates_node(SEXP node) {
    return(VECTOR_ELT(node, SURROGATES_NODE));
}

SEXP info_node(SEXP node) {
    return(VECTOR_ELT(node, INFO_NODE));
}

int is_terminal_node(SEXP node) {
    if (kids_node(node) == R_NilValue) {
        if (split_node(node) != R_NilValue)
            error("kids != split");
        return(1);
    }
    return(0);
}
    
int kidid_node(SEXP node, SEXP data, SEXP vmatch, int obs) {

    SEXP primary, surrogates, prob;
    int ret;
    int i;
    double *dprob;

    primary = split_node(node);
    surrogates = surrogates_node(node);

    /* perform primary split */
    ret = kidid_split(primary, data, vmatch, obs);

    /* surrogate / random splits if needed */
    if (ret == NA_INTEGER) {
        /* surrogate splits */
        if (LENGTH(surrogates) >= 1) {
            for (i = 0; i < LENGTH(surrogates); i++) {
                if (ret != NA_INTEGER) break;
                ret = kidid_split(VECTOR_ELT(surrogates, i), data, vmatch, obs);
            }
        }
        /* random splits */
        if (ret == NA_INTEGER) {
            prob = prob_split(primary);
            dprob = Calloc(LENGTH(prob) - 1, double);
            dprob[0] = REAL(prob)[0];
            for (i = 1; i = LENGTH(prob) - 1; i++)
                dprob[i] = REAL(prob)[i] + dprob[i - 1];
            ret = cut(unif_rand(), dprob, LENGTH(prob) - 1, 1);
            Free(dprob);
        }
    }
    return(ret);
}

int fitted_node(SEXP node, SEXP data, SEXP vmatch, int obs) {

    if (is_terminal_node(node))
        return(id_node(node));
    return(fitted_node(VECTOR_ELT(kids_node(node), 
                                  kidid_node(node, data, vmatch, obs)), 
                      data, vmatch, obs));
}
