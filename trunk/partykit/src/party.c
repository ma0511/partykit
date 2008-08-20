
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

/**
    Linear statistics for conditional inference based on Strasser & Weber (1999)
    *\file LinearStatistic.c
    *\author $Author: hothorn $
    *\date $Date: 2006-08-25 10:53:10 +0200 (Fri, 25 Aug 2006) $
*/
    
void C_abs(double *x, int n) {

    int i;
        for (i = 0; i < n; i++) x[i] = fabs(x[i]);
        }
        
double C_max(const double *x, const int n) {
   double tmp = 0.0;
      int i;
         
            for (i = 0; i < n; i++) {
                   if (x[i] > tmp) tmp = x[i];
                      }
                         return(tmp);
                         }
                         
                         

int NROW (SEXP x) {
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(LENGTH(x));
    return(INTEGER(a)[0]);
}
    
int NCOL (SEXP x) {
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    return(INTEGER(a)[1]);
}

int C_max_factor (SEXP x) {

    int i, max = 0;

    for (i = 0; i < LENGTH(x); i++)
        if (INTEGER(x)[i] > max) max = INTEGER(x)[i];
    return(max);
}

/**
    Computes the Kronecker product of two matrices\n
    *\param A matrix
    *\param m nrow(A)
    *\param n ncol(A)
    *\param B matrix
    *\param r nrow(B)
    *\param s ncol(B)
    *\param ans return value; a pointer to a REALSXP-vector of length (mr x ns)
*/

void C_kronecker (const double *A, const int m, const int n,
                  const double *B, const int r, const int s,
                  double *ans) {

    int i, j, k, l, mr, js, ir;
    double y;

    mr = m * r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++) {
                    ans[(js + l) * mr + ir + k] = y * B[l * r + k];
                }
            }
        }
    }
}  



/**
    Computes the linear statistic, formula (1) in the paper\n
    *\param x values of the transformation
    *\param p dimension of the transformation
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param ans return value; a pointer to a REALSXP-vector of length pq
*/
  
void C_Linstat_numeric (const double *x, const double *y, 
                        const int q, const int *weights, const int n,
                        double *ans) {
              
    int i, k;
    double tmp;

    for (k = 0; k < q; k++)
        ans[k] = 0.0;
            
    for (i = 0; i < n; i++) {
                
        /* optimization: weights are often zero */
        if (weights[i] == 0) continue;
 
        tmp = x[i] * weights[i];
        for (k = 0; k < q; k++)  
            ans[k] += tmp * y[k * n + i];
    }
}


void C_Linstat_factor (const int *x, const int p,
                       const double *y, const int q,
                       const int *weights, const int n,
                       double *ans) {
              
    int i, j, k, kp, kn;
    double tmp;

    for (k = 0; k < q; k++) {

        kn = k * n;
        kp = k * p;
        for (j = 0; j < p; j++) ans[kp + j] = 0.0;
            
        for (i = 0; i < n; i++) {
                
            /* optimization: weights are often zero, treatment contrasts */
            if (weights[i] == 0 || x[i] == 1) continue;
                
            tmp = y[kn + i] * weights[i];
                
            ans[(x[i] - 2) * p + k] += tmp;
        }
    }
}

/**
    Conditional expectation and covariance of the influence function\n
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param ans return value; an object of class `ExpectCovarInfluence'
*/

void C_ExpInf (const double* y, const int q, const int* weights, 
               const int sw, const int n, double* ans) {

    int i, j;    
    
    for (j = 0; j < q; j++) ans[j] = 0.0;
    
    if (sw <= 1) 
        error("C_ExpInf: sum of weights is less than one");

    /*
     *    Expectation of the influence function
     */

    for (i = 0; i < n; i++) {

        /*  observations with zero case weights do not contribute */
    
        if (weights[i] == 0) continue;
    
        for (j = 0; j < q; j++)
            ans[j] += weights[i] * y[j * n + i];
    }

    for (j = 0; j < q; j++)
        ans[j] = ans[j] / sw;

}

/**
    Conditional expectation and covariance of the influence function\n
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param ans return value; an object of class `ExpectCovarInfluence'
*/

void C_CovInf (const double* y, const int q, const int* weights, const int sw, 
               const int n, const double* exp, double *ans) {

    int i, j, k, jq;
    double tmp;
    
    for (j = 0; j < q*q; j++) ans[j] = 0.0;

    if (sw <= 1) 
        error("C_CovInf: sum of weights is less than one");
    
    /*
     *    Covariance of the influence function
     */

    for (i = 0; i < n; i++) {

        if (weights[i] == 0) continue;
     
        for (j = 0; j < q; j++) {
            tmp = weights[i] * (y[j * n + i] - exp[j]);
            jq = j * q;
            for (k = 0; k < q; k++)
                ans[jq + k] += tmp * (y[k * n + i] - exp[k]);
        }
    }

    for (j = 0; j < q*q; j++)
        ans[j] = ans[j] / sw;
}

void C_swx_numeric (const double *x, const int *weights, const int n, 
                    double *swx, double *swx2) {

    int i;
    double tmp;

    swx[0] = 0.0;
    swx2[0] = 0.0;

    for (i = 0; i < n; i++) {

        /*  observations with zero case weights do not contribute */
        if (weights[i] == 0) continue;
    
        tmp = weights[i] * x[i];
        swx[0] += tmp;

        /* covariance part */
        swx2[0] += tmp * x[i];
    }
}

void C_swx_factor (const int *x, const int p, const int *weights, const int n,
                   double *swx, double *swx2) {

    int i, j, k;

    for (j = 0; j < p; j++) {
        swx[j] = 0.0;
        for (k = 0; k < p; k++)
            swx2[j * p + k] = 0.0;
    }

    for (i = 0; i < n; i++) {

        /*  observations with zero case weights do not contribute, treatment contrasts */
        if (weights[i] == 0 || x[i] == 1) continue;
    
        swx[x[i] - 2] += (double) weights[i];
    }
    
    /* covariance part */
    for (j = 0; j < p; j++)
        swx[j * p + j] = swx[j];
}


/**
    Conditional expectation and covariance of the a linear statistic\n
    *\param x values of the transformation
    *\param p dimension of the transformation
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param expcovinf an object of class `ExpectCovarInfluence'
    *\param ans return value; an object of class `ExpectCovar'
*/

void C_ExpCovLinstat (const double *swx, const double *swx2, const int p, const int q, 
                      const int sw, const double *expinf, const double *covinf, 
                      double *explinstat, double *covlinstat) {

    int i, j, k, pq, ip;
    double f1, f2, tmp, dsw;
    double *CT2, *Covy_x_swx;

    pq = p * q;

    /*
    *   explinstat: expectation of the linear statistic T
    */

    for (k = 0; k < p; k++) {
        for (j = 0; j < q; j++)
            explinstat[j * p + k] = swx[k] * expinf[j];
    }

    if (sw <= 1.0) 
        error("C_ExpCovLinstat: sum of weights is less than one");

    /* 
    *   covlinstat:  covariance of the linear statistic T
    */

    dsw = (double) sw;
    f1 = dsw / (dsw - 1);
    f2 = (1 / (dsw - 1));

    if (pq == 1) {
        covlinstat[0] = f1 * covinf[0] * swx2[0];
        covlinstat[0] -= f2 * covinf[0] * swx[0] * swx[0];
    } else {
        /* two more helpers needed */
        CT2 = Calloc(pq * pq, double);            /* pq x pq */
        Covy_x_swx = Calloc(pq * q, double);      /* pq x q  */
        
        C_kronecker(covinf, q, q, swx2, p, p, covlinstat);
        C_kronecker(covinf, q, q, swx, p, 1, Covy_x_swx);
        C_kronecker(Covy_x_swx, pq, q, swx, 1, p, CT2);

        for (k = 0; k < (pq * pq); k++)
            covlinstat[k] = f1 * covlinstat[k] - f2 * CT2[k];

        /* clean up */
        Free(CT2); Free(Covy_x_swx);
    }
}


/**
    R-interface to C_ExpCovLinstat\n
    *\param x values of the transformation
    *\param y values of the influence function
    *\param weights case weights
    *\param expcovinf an object of class `ExpectCovarInfluence'
*/

SEXP C_LinstatExpCov(SEXP x, SEXP y, SEXP weights, SEXP ans) {
    
    SEXP expcovinf, explinstat, covlinstat, linstat, dim, dx;
    double *expinf, *covinf, *swx, *swx2;
    int sw = 0, i;
    int n, p, q, pq;

    /* determine the dimensions and some checks */

    n  = NROW(x);
    q  = NCOL(y);
    
    if (NROW(y) != n)
        error("R_ExpCovLinstat: y does not have %d rows", n);
    if (LENGTH(weights) != n) 
        error("R_ExpCovLinstat: vector of case weights does not have %d elements", n);

    for (i = 0; i < n; i++) sw += INTEGER(weights)[i];

    if (isReal(x)) {
        p = 1;
        swx = Calloc(p, double);
        swx2 = Calloc(p, double);
        
        C_swx_numeric(REAL(x), INTEGER(weights), n, swx, swx2);
    } 
    if (isFactor(x) & !isOrdered(x)) {
        p = C_max_factor(x) - 1;
        swx = Calloc(p, double);
        swx2 = Calloc(p, double);
                        
        C_swx_factor(INTEGER(x), p, INTEGER(weights), n, swx, swx2);                      
    }
    /* FIXME: implement scores for ordered */
    if (isInteger(x) || isOrdered(x)) {
        PROTECT(dx = coerceVector(x, REALSXP));
        p = 1;
        swx = Calloc(p, double);
        swx2 = Calloc(p, double);
                        
        C_swx_numeric(REAL(dx), INTEGER(weights), n, swx, swx2);
        UNPROTECT(1);
    }
                                
        
    pq = p * q;

    expinf = Calloc(q, double);
    covinf = Calloc(q * q, double);

    C_ExpInf(REAL(y), q, INTEGER(weights), sw, n, expinf);
    C_CovInf(REAL(y), q, INTEGER(weights), sw, n, expinf,  covinf);

    SET_VECTOR_ELT(ans, 0, dim = allocVector(INTSXP, 2));    
    SET_VECTOR_ELT(ans, 1, linstat = allocVector(REALSXP, pq));
    SET_VECTOR_ELT(ans, 2, explinstat = allocVector(REALSXP, pq));
    SET_VECTOR_ELT(ans, 3, covlinstat = allocVector(REALSXP, pq * pq));
    INTEGER(dim)[0] = p;
    INTEGER(dim)[1] = q;
     
    if (isReal(x)) {
        C_Linstat_numeric(REAL(x), REAL(y), q, INTEGER(weights), n, REAL(linstat));
    } else {
        C_Linstat_factor(INTEGER(x), p, REAL(y), q, INTEGER(weights), n, REAL(linstat));
    }

    C_ExpCovLinstat(swx, swx2, p, q, sw, expinf, covinf, REAL(explinstat), 
                    REAL(covlinstat));
    
    return(ans);
}

SEXP R_LinstatExpCov(SEXP x, SEXP y, SEXP weights) {

    SEXP ans;
    
    PROTECT(ans = allocVector(VECSXP, 4));
    C_LinstatExpCov(x, y, weights, ans);
    UNPROTECT(1);
    return(ans);
}


/**
    Standardizes a statistic t of length pq with mean mu and covariance Sigma
    for variances > tol \n
    *\param t the vector of statistics
    *\param mu expectations
    *\param Sigma covariance matrix
    *\param pq dimension of t
    *\param tol tolerance for variances
    *\param ans return value; a pointer to a REALSXP-vector of length pq
*/

void C_standardize (const double *t, const double *mu, const double *Sigma, 
                    int pq, double tol, double *ans) {
                   
    int i;
    double sd;
    
    for (i = 0; i < pq; i++) {
        sd = Sigma[i*pq + i]; 
        if (sd > tol)
            ans[i] = (t[i] - mu[i])/sqrt(sd);
        else
            ans[i] = 0.0;
    }
}


/**
    Absolute values of standardized statistics
    *\param t the vector of statistics
    *\param mu expectations
    *\param Sigma covariance matrix
    *\param pq dimension of t
    *\param tol tolerance for variances
    *\param ans return value; a pointer to a REALSXP-vector of length pq
*/

void C_absstandardize (const double *t, const double *mu, const double *Sigma,
                       int pq, double tol, double *ans) {
                      
    C_standardize(t, mu, Sigma, pq, tol, ans);
    C_abs(ans, pq);
}


/**
    Maximum absolute values of standardized statistics
    *\param t the vector of statistics
    *\param mu expectations
    *\param Sigma covariance matrix
    *\param pq dimension of t
    *\param tol tolerance for variances
*/

double C_maxabsTestStatistic (const double *t, const double *mu, const double *Sigma,
                              int pq, double tol) {
                           
     double *mem, ans;
     
     mem = Calloc(pq, double);
     C_absstandardize(t, mu, Sigma, pq, tol, mem);
     ans = C_max(mem, pq);
     Free(mem);
     return(ans);
}

double C_split_numeric (SEXP x, SEXP y, SEXP weights, int minsplit) {

     int *xsort, i, j, k, n, sw, q, *iweights;
     double dsw = 0.0, *expinf, *covinf, *dx, tmp, lastx, thisx;
     double breaks = 0.0, cstat = 0.0, laststat = 0.0;
     double *linstat, *explinstat, *covlinstat, *dy;

     n = NROW(x);
     if (NCOL(x) != 1 | !isReal(x))
         error("");
         
     xsort = Calloc(n, int);
     for (i = 0; i < n; i++) xsort[i] = i;

     dx = Calloc(n, double);
     for (i = 0; i < n; i++) dx[i] = REAL(x)[i];
     rsort_with_index(dx, xsort, n);
     Free(dx);
     dx = REAL(x);
    
     iweights = INTEGER(weights);
     for (i = 0; i < n; i++) sw += iweights[i];

     q = NCOL(y);     
     dy = REAL(y);
     expinf = Calloc(q, double);
     covinf = Calloc(q * q, double);

     C_ExpInf(REAL(y), q, iweights, sw, n, expinf);
     C_CovInf(REAL(y), q, iweights, sw, n, expinf,  covinf);
     
     linstat = Calloc(q, double);
     explinstat = Calloc(q, double);
     covlinstat = Calloc(q * q, double);
     
     for (j = 0; j < q; j++) {
         linstat[j] = 0.0;
         explinstat[j] = 0.0;
         for (k = 0; k < q; k++) 
             covlinstat[j * q + k] = 0.0;
     }
     
     dsw = 0;
     lastx = NA_REAL;
     for (i = 0; i < n; i++) {

         thisx = dx[xsort[i]];
         if (iweights[xsort[i]] == 0 | thisx == NA_REAL) continue;
         
         dsw += (double) iweights[xsort[i]];
         for (k = 0; k < q; k++)
             linstat[k] += iweights[xsort[i]] * dy[k * n + xsort[i]];

         if (dsw > (n - minsplit)) break;
         if (dsw < minsplit) continue;
         if (thisx == lastx) continue;
         lastx = thisx;
             
         C_ExpCovLinstat(&dsw, &dsw, 1, q, sw, expinf, covinf, explinstat,
                         covlinstat);

         cstat = C_maxabsTestStatistic(linstat, explinstat, covlinstat, q, 0.0);

         if (cstat > laststat) {
             breaks = lastx;
             laststat = cstat;
         }
     }

     Free(expinf);
     Free(covinf);
     Free(linstat);
     Free(explinstat);
     Free(covlinstat);

     return(breaks);
}

SEXP R_split_numeric (SEXP x, SEXP y, SEXP weights, SEXP minsplit) {

    SEXP ans, dx;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    
    if (isInteger(x) | isOrdered(x)) {
        PROTECT(dx = coerceVector(x, REALSXP));
        REAL(ans)[0] = C_split_numeric(dx, y, weights, INTEGER(minsplit)[0]);
        UNPROTECT(1);
    } else {
        if (isFactor(x)) error("x is a factor");
        REAL(ans)[0] = C_split_numeric(x, y, weights, INTEGER(minsplit)[0]);
    }

    UNPROTECT(1);
    return(ans);
}

SEXP R_split_factor (SEXP x, SEXP y, SEXP weights, int minsplit) {

    SEXP tmpx, lec, ans;
    int i, jj, mi, l, j, p, q, *indl, n;
    int *iweights, *ix, sw;
    double *expinf, *covinf,  *dtmpx, *linstat, *explinstat, *covlinstat;
    double cstat = 0.0, laststat = 0.0;
    

    n = NROW(x);
    PROTECT(tmpx = allocVector(REALSXP, n));
    dtmpx = REAL(tmpx);
    ix = INTEGER(x);

    PROTECT(lec = allocVector(VECSXP, 4));

    iweights = INTEGER(weights);
    for (i = 0; i < n; i++) sw += iweights[i];

    q = NCOL(y);     
    expinf = Calloc(q, double);
    covinf = Calloc(q * q, double);

    C_ExpInf(REAL(y), q, iweights, sw, n, expinf);
    C_CovInf(REAL(y), q, iweights, sw, n, expinf, covinf);
     
    
    laststat = 0.0;
    sw = 0;
    mi = 1;  
    p = C_max_factor(x);
    PROTECT(ans = allocVector(INTSXP, p));
    indl = Calloc(p, int);

    for (l = 1; l < p; l++) mi *= 2;
    for (j = 0; j < mi-1; j++) { /* go though all splits */
        jj = j;
        for (l = 1; l < p; l++) {
            indl[l] = (jj%2);
            jj /= 2;
        }
        
        for (i = 0; i < n; i++) dtmpx[i] = 0.0;
        sw = 0;
        
        for (i = 0; i < n; i++) {
        
            if (iweights[i] > 0 & indl[ix[i] - 1]) {
                dtmpx[i] += 1.0;
                sw++;
            }
        }
        if (sw > (n - minsplit)) continue;
        if (sw < minsplit) continue;

        C_LinstatExpCov(tmpx, y, weights, lec);
        linstat = REAL(VECTOR_ELT(lec, 1));
        explinstat = REAL(VECTOR_ELT(lec, 2));
        covlinstat = REAL(VECTOR_ELT(lec, 3));
        
        cstat = C_maxabsTestStatistic(linstat, explinstat, covlinstat, q, 0.0);
        if (cstat > laststat) {
            laststat = cstat;
            for (l = 1; l < p; l++)
                INTEGER(ans)[l] = indl[l] + 1;
        }
    }
  
    UNPROTECT(3);
    return(ans);
}                                          
