
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>


#define VARID_SPLIT		0		/* varid_split */
#define BREAKS_SPLIT		1		/* breaks_split */
#define INDEX_SPLIT		2		/* index_split */
#define RIGHT_SPLIT		3		/* right_split */
#define PROB_SPLIT		4		/* prob_split */
#define INFO_SPLIT		5		/* info_split */
#define LENGTH_SPLIT		INFO_SPLIT + 1	/* length split object */

void C_split(SEXP varid, SEXP breaks, SEXP index, SEXP right,
             SEXP prob, SEXP info, SEXP split);

int varid_split(SEXP split);
SEXP breaks_split(SEXP split);
SEXP index_split(SEXP split);
int right_split(SEXP split);
SEXP prob_split(SEXP split);
SEXP info_split(SEXP split);
