
#include "partykit.h"

int cut(double x, double *breaks, int n, int right) {

    int ret, i;

    ret = NA_INTEGER;
    if (ISNA(x)) return(ret);

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
