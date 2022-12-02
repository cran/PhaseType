#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "utility.h"


double *LJMA_LAPACK_work;
int LJMA_LAPACK_lwork = -1;


int LJMA_counter = 0;


// The next two functions are for allocating the LAPACK work space when calling certain functions directly from R instead of via the outer Gibbs call
void LJMA_LAPACKspace(int *n) {
	char balanc = 'B', jobvl = 'V', jobvr = 'V', sense = 'B'; int lwork = -1, info; double work;
	F77_CALL(dgeevx)(&balanc, &jobvl, &jobvr, &sense, n, NULL, n, NULL, NULL, NULL, n, NULL, n, NULL, NULL, NULL, NULL, NULL, NULL, &work, &lwork, NULL, &info FCONE FCONE FCONE FCONE);
	LJMA_LAPACK_lwork = (int) work;
	F77_CALL(dgetri)(n, NULL, n, NULL, &work, &lwork, &info);
	if((int) work > LJMA_LAPACK_lwork) LJMA_LAPACK_lwork = (int) work;
	//LJMA_LAPACK_work = (double *) R_alloc(LJMA_LAPACK_lwork, sizeof(double));
	LJMA_LAPACK_work = (double *) Calloc(LJMA_LAPACK_lwork, double);
}
void LJMA_LAPACKspaceFree(void) {
	Free(LJMA_LAPACK_work);
}


// Then we can get at the performance counter here
void LJMA_setCounter(int *counter) {
	LJMA_counter = *counter;
}
void LJMA_getCounter(int *counter) {
	*counter = LJMA_counter;
}


// Functions to track random number generator state
int LJMA_RNG = 0;


//// Invert a square matrix
// A (input/output)
//     square matrix to be inverted ... overwritten with the inverse
// n (input)
//     dimension
// workI (output)
//     n element workspace
int LJMA_inverse(double *A, int *n, int *workI) {
	// Do LU decomposition as required by LAPACK for inverse calculation
	int *ipiv, info;
	ipiv = workI;
	F77_CALL(dgetrf)(n, n, A, n, ipiv, &info);
		if(info != 0) {
			Rprintf("Error (LJMA_inverse 01): failed LAPACK call, code=%d\n", info);
			return(info);
		}

	// Do inverse calculation using that LU decomposition
	F77_CALL(dgetri)(n, A, n, ipiv, LJMA_LAPACK_work, &LJMA_LAPACK_lwork, &info);
		if(info != 0) {
			Rprintf("Error (LJMA_inverse 03): failed LAPACK call, code=%d\n", info);
			return(info);
		}
	return(0);
}


//// Compute the eigen decomposition of a matrix
// n (input)
//     dimension of the matrix
// S (input)
//     matrix to decompose
// evals (output)
//     array of length n in which to store eigenvalues
// Q (output)
//     n x n matrix in which to store eigenvectors
// Qinv (output)
//     n x n matrix in which to store the inverse of the eigenvectors matrix
// workD (output)
//     2n^2 + 4n element workspace
// workI (output)
//     3n - 2 element workspace
//
// return: 0 on success or non-0 on failure
int LJMA_eigen(int *n, double *S, double *evals, double *Q, double *Qinv, double *workD, int *workI) {
	// copy S as we want to be non-destructive, but LAPACK routines will muller it
	double *A;
	A = workD; workD += *n * *n;
	memcpy(A, S, *n * *n * sizeof(*S));

	// Setup storage for imaginary part of eigenvalues; the left eigenvectors; the scaling details (even though unwanted) and condition numbers and workspace
	double *evalsi;
	evalsi = workD; workD += *n;
	double *Ql;
	Ql = workD; workD += *n * *n;
	double *scale;
	scale = workD; workD += *n;
	double *rconde; // reciprocal of e-val condition numbers
	rconde = workD; workD += *n;
	double *rcondv; // reciprocal of e-vect condition numbers
	rcondv = workD; workD += *n;
	int *iwork;
	iwork = workI; workI += 2 * *n - 2;

	// do eigen decomposition
	char balanc = 'B', jobvl = 'V', jobvr = 'V', sense = 'B'; int info, ilo, ihi; double abnrm;
	F77_CALL(dgeevx)(&balanc, &jobvl, &jobvr, &sense, n, A, n, evals, evalsi, Ql, n, Q, n, &ilo, &ihi, scale, &abnrm, rconde, rcondv, LJMA_LAPACK_work, &LJMA_LAPACK_lwork, iwork, &info FCONE FCONE FCONE FCONE);
		if(info != 0) {
			Rprintf("Error (LJMA_eigen 01): failed LAPACK call, code=%d\n", info);
			return(info);
		}

	// Check none of the eigen values have an imaginary part & condition numbers are ok
	//char cmach = 'E';
	//double precision = F77_CALL(dlamch)(&cmach FCONE);
	for(int i=0; i<*n; i++) {
		if(evalsi[i] > 0) Rprintf("Error: imaginary part of eigenvalue %d found.\n", i+1);
		//Rprintf("i=%d: eval = %e, err(eval) = %e, err(evect) = %e\n", i+1, evals[i], precision*abnrm/rconde[i], precision*abnrm/rcondv[i]);
	}

	//// Now find Q inverse
	// copy Q into Qinv
	memcpy(Qinv, Q, *n * *n * sizeof(*Q));
	LJMA_inverse(Qinv, n, workI);

	return(0);
}




/* The following was taken from core R for personal
 * maintenance/modification.
 * It is was private non-API code as noted by Prof Ripley in personal
 * correspondence, so calls to it were removed from this package and
 * the source brought into the package for use directly.
 */
/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999, 2001 the R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* from NETLIB c/brent.shar with max.iter, add'l info and convergence
 details hacked in by Peter Dalgaard */

/*************************************************************************
 *			    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,info,Tol,Maxit)
 *	double ax;			Root will be seeked for within
 *	double bx;			a range [ax,bx]
 *	double (*f)(double x, void *info); Name of the function whose zero
 *					will be seeked for
 *	void *info;			Add'l info passed to f
 *	double *Tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 *	int *Maxit;			Max. iterations
 *
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *	*Tol returns estimated precision
 *	*Maxit returns actual # of iterations, or -1 if maxit was
 *      reached without convergence.
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bisection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *		   the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bisection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bisection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */

#define EPSILON DBL_EPSILON

double Find0(			/* An estimate of the root */
			 double ax,				/* Left border | of the range	*/
			 double bx,				/* Right border| the root is seeked*/
			 double (*f)(double x, void *info),	/* Function under investigation	*/
			 void *info,				/* Add'l info passed on to f	*/
			 double *Tol,			/* Acceptable tolerance		*/
			 int *Maxit)				/* Max # of iterations */
{
    double fa = (*f)(ax, info);
    double fb = (*f)(bx, info); //Rprintf("L(%lf): %e, R(%lf): %e\n",ax,fa,bx,fb);
    return Find02(ax, bx, fa, fb, f, info, Tol, Maxit);
}

/* faster for "expensive" f(), in those typical cases where
 *             f(ax) and f(bx) are available anyway : */

double Find02(			/* An estimate of the root */
			  double ax,				/* Left border | of the range	*/
			  double bx,				/* Right border| the root is seeked*/
			  double fa, double fb,		/* f(a), f(b) */
			  double (*f)(double x, void *info),	/* Function under investigation	*/
			  void *info,				/* Add'l info passed on to f	*/
			  double *Tol,			/* Acceptable tolerance		*/
			  int *Maxit)				/* Max # of iterations */
{
    double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
    double tol;
    int maxit;

    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = * Tol;

    /* First test if we have found a root at an endpoint */
    if(fa == 0.0) {
		*Tol = 0.0;
		*Maxit = 0;
		return a;
    }
    if(fb ==  0.0) {
		*Tol = 0.0;
		*Maxit = 0;
		return b;
    }

    while(maxit--)		/* Main iteration loop	*/
    {
		double prev_step = b-a;		/* Distance from the last but one
									 to the last approximation	*/
		double tol_act;			/* Actual tolerance		*/
		double p;			/* Interpolation step is calcu- */
		double q;			/* lated in the form p/q; divi-
							 * sion operations is delayed
							 * until the last moment	*/
		double new_step;		/* Step at this iteration	*/

		if( fabs(fc) < fabs(fb) )
		{				/* Swap data for b to be the	*/
			a = b;  b = c;  c = a;	/* best approximation		*/
			fa=fb;  fb=fc;  fc=fa;
		}
		tol_act = 2*EPSILON*fabs(b) + tol/2;
		new_step = (c-b)/2;

		if( fabs(new_step) <= tol_act || fb == (double)0 )
		{
			*Maxit -= maxit;
			*Tol = fabs(c-b);
			return b;			/* Acceptable approx. is found	*/
		}

		/* Decide if the interpolation can be tried	*/
		if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
		   && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
										 * Interpolation may be tried	*/
			register double t1,cb,t2;
			cb = c-b;
			if( a==c ) {		/* If we have only two distinct	*/
				/* points linear interpolation	*/
				t1 = fb/fa;		/* can only be applied		*/
				p = cb*t1;
				q = 1.0 - t1;
			}
			else {			/* Quadric inverse interpolation*/

				q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
				p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);
			}
			if( p>(double)0 )		/* p was calculated with the */
				q = -q;			/* opposite sign; make p positive */
			else			/* and assign possible minus to	*/
				p = -p;			/* q				*/

			if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
			   && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
				new_step = p/q;			/* it is accepted
										 * If p/q is too large then the
										 * bisection procedure can
										 * reduce [b,c] range to more
										 * extent */
		}

		if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
			if( new_step > (double)0 )	/* than tolerance		*/
				new_step = tol_act;
			else
				new_step = -tol_act;
		}
		a = b;	fa = fb;			/* Save the previous approx. */
		b += new_step;	fb = (*f)(b, info);	/* Do step to a new approxim. */
		if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
			/* Adjust c for it to have a sign opposite to that of b */
			c = a;  fc = fa;
		}

    }
    /* failed! */
    *Tol = fabs(c-b);
    *Maxit = -1;
    return b;
}
