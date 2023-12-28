#include "f2c.h"
#include "blaswrap.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int sorg2l_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    integer i__, j, l, ii;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    slarf_(char *, integer *, integer *, real *, integer *, real *, 
	    real *, integer *, real *), xerbla_(char *, integer *);


/*  -- LAPACK routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SORG2L generates an m by n real matrix Q with orthonormal columns, */
/*  which is defined as the last n columns of a product of k elementary */
/*  reflectors of order m */

/*        Q  =  H(k) . . . H(2) H(1) */

/*  as returned by SGEQLF. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix Q. M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix Q. M >= N >= 0. */

/*  K       (input) INTEGER */
/*          The number of elementary reflectors whose product defines the */
/*          matrix Q. N >= K >= 0. */

/*  A       (input/output) REAL array, dimension (LDA,N) */
/*          On entry, the (n-k+i)-th column must contain the vector which */
/*          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/*          returned by SGEQLF in the last k columns of its array */
/*          argument A. */
/*          On exit, the m by n matrix Q. */

/*  LDA     (input) INTEGER */
/*          The first dimension of the array A. LDA >= max(1,M). */

/*  TAU     (input) REAL array, dimension (K) */
/*          TAU(i) must contain the scalar factor of the elementary */
/*          reflector H(i), as returned by SGEQLF. */

/*  WORK    (workspace) REAL array, dimension (N) */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit */
/*          < 0: if INFO = -i, the i-th argument has an illegal value */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SORG2L", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

/*     Initialise columns 1:n-k to columns of the unit matrix */

    i__1 = *n - *k;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    a[l + j * a_dim1] = 0.f;
/* L10: */
	}
	a[*m - *n + j + j * a_dim1] = 1.f;
/* L20: */
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *n - *k + i__;

/*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left */

	a[*m - *n + ii + ii * a_dim1] = 1.f;
	i__2 = *m - *n + ii;
	i__3 = ii - 1;
	slarf_("Left", &i__2, &i__3, &a[ii * a_dim1 + 1], &c__1, &tau[i__], &
		a[a_offset], lda, &work[1]);
	i__2 = *m - *n + ii - 1;
	r__1 = -tau[i__];
	sscal_(&i__2, &r__1, &a[ii * a_dim1 + 1], &c__1);
	a[*m - *n + ii + ii * a_dim1] = 1.f - tau[i__];

/*        Set A(m-k+i+1:m,n-k+i) to zero */

	i__2 = *m;
	for (l = *m - *n + ii + 1; l <= i__2; ++l) {
	    a[l + ii * a_dim1] = 0.f;
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     End of SORG2L */

} /* sorg2l_ */
