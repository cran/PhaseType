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
	F77_CALL(dgeevx)(&balanc, &jobvl, &jobvr, &sense, n, NULL, n, NULL, NULL, NULL, n, NULL, n, NULL, NULL, NULL, NULL, NULL, NULL, &work, &lwork, NULL, &info);
	LJMA_LAPACK_lwork = (int) work;
	F77_CALL(dgetri)(n, NULL, n, NULL, &work, &lwork, &info);
	if((int) work > LJMA_LAPACK_lwork) LJMA_LAPACK_lwork = (int) work;
	//LJMA_LAPACK_work = (double *) R_alloc(LJMA_LAPACK_lwork, sizeof(double));
	LJMA_LAPACK_work = (double *) Calloc(LJMA_LAPACK_lwork, double);
}
void LJMA_LAPACKspaceFree() {
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
	F77_CALL(dgeevx)(&balanc, &jobvl, &jobvr, &sense, n, A, n, evals, evalsi, Ql, n, Q, n, &ilo, &ihi, scale, &abnrm, rconde, rcondv, LJMA_LAPACK_work, &LJMA_LAPACK_lwork, iwork, &info);
		if(info != 0) {
			Rprintf("Error (LJMA_eigen 01): failed LAPACK call, code=%d\n", info);
			return(info);
		}
	
	// Check none of the eigen values have an imaginary part & condition numbers are ok
	//char cmach = 'E';
	//double precision = F77_CALL(dlamch)(&cmach);
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
