#include "f2c.h"
#include "blaswrap.h"

/* Table of constant values */

static integer c__7 = 7;
static doublereal c_b5 = 0.;
static doublereal c_b6 = 1.;

doublereal drzt02_(integer *m, integer *n, doublereal *af, integer *lda, 
	doublereal *tau, doublereal *work, integer *lwork)
{
    /* System generated locals */
    integer af_dim1, af_offset, i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    integer i__, info;
    doublereal rwork[1];
    extern doublereal dlamch_(char *), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *), 
	    xerbla_(char *, integer *), dormrz_(char *, char *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *, 
	     integer *);


/*  -- LAPACK test routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DRZT02 returns */
/*       || I - Q'*Q || / ( M * eps) */
/*  where the matrix Q is defined by the Householder transformations */
/*  generated by DTZRZF. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix AF. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix AF. */

/*  AF      (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*          The output of DTZRZF. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array AF. */

/*  TAU     (input) DOUBLE PRECISION array, dimension (M) */
/*          Details of the Householder transformations as returned by */
/*          DTZRZF. */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK) */

/*  LWORK   (input) INTEGER */
/*          length of WORK array. LWORK >= N*N+N*NB. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    af_dim1 = *lda;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    --tau;
    --work;

    /* Function Body */
    ret_val = 0.;

    if (*lwork < *n * *n + *n) {
	xerbla_("DRZT02", &c__7);
	return ret_val;
    }

/*     Quick return if possible */

    if (*m <= 0 || *n <= 0) {
	return ret_val;
    }

/*     Q := I */

    dlaset_("Full", n, n, &c_b5, &c_b6, &work[1], n);

/*     Q := P(1) * ... * P(m) * Q */

    i__1 = *n - *m;
    i__2 = *lwork - *n * *n;
    dormrz_("Left", "No transpose", n, n, m, &i__1, &af[af_offset], lda, &tau[
	    1], &work[1], n, &work[*n * *n + 1], &i__2, &info);

/*     Q := P(m) * ... * P(1) * Q */

    i__1 = *n - *m;
    i__2 = *lwork - *n * *n;
    dormrz_("Left", "Transpose", n, n, m, &i__1, &af[af_offset], lda, &tau[1], 
	     &work[1], n, &work[*n * *n + 1], &i__2, &info);

/*     Q := Q - I */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[(i__ - 1) * *n + i__] += -1.;
/* L10: */
    }

    ret_val = dlange_("One-norm", n, n, &work[1], n, rwork) / (
	    dlamch_("Epsilon") * (doublereal) max(*m,*n));
    return ret_val;

/*     End of DRZT02 */

} /* drzt02_ */
