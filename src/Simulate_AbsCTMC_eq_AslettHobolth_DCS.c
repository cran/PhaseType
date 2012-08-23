#include <R.h>
#include <Rmath.h>

#include <R_ext/BLAS.h>
#include "utility.h"
#include "Simulate_AbsCTMC_eq_AslettHobolth_DCS.h"
#include "Simulate_AbsCTMC_gt_Hobolth_DCS.h"


/* pi %*% Q %*% exp(L*y) %*% Qinv * s_i for each i */
int LJMA_Hobolth_endState(double y, double *pi, double *Q, double *evals, double *Qinv, double *s, int n, double *workD) {
	LJMA_GetRNGstate();
	const char transT = 'T'; const double oneD = 1.0, zeroD = 0.0; const int oneI = 1.0;
	
	double *p, *tmp;
	p = workD; workD += n;
	tmp = workD; workD += n;
	
	// p <- pi %*% Q
	F77_CALL(dgemv)(&transT, &n, &n, &oneD, Q, &n, pi, &oneI, &zeroD, p, &oneI);
	// p <- pi %*% Q %*% exp(L*y)
	for(int i=0; i<n; i++) {
		p[i] *= exp(evals[i]*y);
	}
	// tmp <- pi %*% Q %*% exp(L*y) %*% Qinv
	F77_CALL(dgemv)(&transT, &n, &n, &oneD, Qinv, &n, p, &oneI, &zeroD, tmp, &oneI);
	// p <- pi %*% Q %*% exp(L*y) %*% Qinv * s
	double sum = 0.0;
	for(int i=0; i<n; i++) {
		sum += p[i] = tmp[i] * s[i];
	}
	for(int i=0; i<n; i++) {
		p[i] = p[i]/sum;
		//Rprintf("Prob -> %d = %lf\n", i+1, p[i]);
	}
	if(p[0]!=0.0) {
		Rprintf("ALERT!\n");
		Rprintf("y=%lf\npi=(%lf,%lf,%lf)\nQ={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}}\n", y, pi[0], pi[1], pi[2], Q[0], Q[1], Q[2], Q[3], Q[4], Q[5], Q[6], Q[7], Q[8]);
	}
	
	double sofar = 0.0, target;
	target = runif(0.0, 1.0);
	int res = 0;
	while(sofar < target) {
		sofar += p[res++];
	}
	res--;
	
	LJMA_PutRNGstate();
	return(res);
}


//// Aslett modified Hobolth sampling of an absorbing chain using directly conditional sampling, returning only sufficient statistics of the chain
// y (input)
//     m dimensional array of absorption times to be conditioned on
// censored (input)
//     m dimensional array of 1/0 for is/not a censored observation
// m (input)
//     number of observations
// pi (input)
//     1 x n vector of initial state probabilities (*must* sum to 1)
// S (input)
//     n x n matrix of transition rates between non-absorbing states
// s (input)
//     1 x n vector of exit rates
// Q (input)
//     precomputed n x n matrix of eigenvectors of the original S matrix
// evals (input)
//     array size n of the n eigenvalues of the original S matrix
// Qinv_b (input)
//     precomputed n x 1 matrix Qinv %*% b, where b is vector with 1 for each exitable state
// b (input)
//     1 x n vector of 0/1s where 1 indicates an exit to be allowed
// Qinv (input)
//     precomputed n x n matrix Q^{-1} [NB ignored unless reverse=1]
// n (input)
//     dimension
// iter (input)
//     number of MH iterations to perform
// res_z (output)
//     n dimensional array, filled with sums of times spent in each state
// res_B (output)
//     n dimensional array, set to the number of times chains started in the given states
// res_N (output)
//     n x n matrix, filled with the sums of the number of transitions between states (i != j).  Diagonal entries are the exit-to-absorption counts
// workD (output)
//     ?? element workspace
// workI (output)
//     ?? element workspace
//
void LJMA_MHsample_Hobolth2(double *y, int *censored, int *m, double *pi, double *S, double *s, double *Q, double *evals, double *Qinv_b, double *bvec, double *Qinv, int *n, int *iter, double *res_z, int *res_B, int *res_N, double *workD, int *workI) {
	// Initialise output vars
	int i,j,k,l;
	for(i=0; i<*n; i++) {
		res_B[i] = 0;
		res_z[i] = 0.0;
		for(j=0; j<*n; j++) {
			res_N[i + j * *n] = 0;
		}
	}
	
	// Create memory for current steps
	double *c_z;
	int c_B, *c_N, c_pre;
	c_z = workD; workD += *n;
	c_N = workI; workI += *n * *n;
	
	// Exit state choice
	int b; //double *bvec;
	//bvec = workD; workD += *n;
	const char transN = 'N'; const double oneD = 1.0, zeroD = 0.0; const int oneI = 1.0;
	
	// Off we go
	double *y_p;
	int *censored_p;
	y_p = y;
	censored_p = censored;
	for(i=0; i<*m; i++) {
		// Decide what state to end in before running sample ...
		b = LJMA_Hobolth_endState(*y_p, pi, Q, evals, Qinv, s, *n, workD);
			//Rprintf("b=%d (%lf,%lf,%lf)\n", b, s[0],s[1],s[2]);
		// ... and setup Qinv_b
		for(j=0; j<*n; j++) {
			if(j == b)
				bvec[j] = 1.0;
			else
				bvec[j] = 0.0;
		}
		F77_CALL(dgemv)(&transN, n, n, &oneD, Qinv, n, bvec, &oneI, &zeroD, Qinv_b, &oneI);
		
		if(*censored_p) Rprintf("\nWARNING: DCS does not support censoring\n");
		LJMA_samplechain_Hobolth(y_p, pi, S, Q, evals, Qinv_b, bvec, s, n, c_z, &c_B, c_N, &c_pre, workD);
		
		// Now add the current step to the running total
		res_B[c_B]++;
		for(k=0; k<*n; k++) {
			res_z[k] += c_z[k];
			for(l=0; l<*n; l++) {
				res_N[k + l * *n] += c_N[k + l * *n];
			}
		}
		
		y_p++;
		censored_p++;
	}
}

