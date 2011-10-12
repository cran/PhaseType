#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "Simulate_AbsCTMC_eq_Aslett_DCS.h"
#include "Simulate_AbsCTMC_gt_Aslett_DCS.h"


//// (ASLETT METHOD) Metropolis-Hastings sample a chain conditional on absorbing after an observation, returning only sufficient statistics of the chain
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
// P (input)
//     n x n matrix of embedded transition probabilities, excluding absorbing moves
// Pfull (input)
//     n x (n+1) matrix of embedded transition probabilities
// Q (input)
//     precomputed n x n matrix of eigenvectors of the original S matrix
// evals (input)
//     array size n of the n eigenvalues of the original S matrix
// Qinv_1 (input)
//     precomputed n x 1 matrix Qinv %*% 1
// Qinv (input)
//     precomputed n x n matrix Q^{-1} [NB ignored unless reverse=1]
// n (input)
//     dimension
// iter (input)
//     number of MH iterations to perform
// reverse (input)
//     sample in reverse
// res_z (output)
//     n dimensional array, filled with sums of times spent in each state
// res_B (output)
//     n dimensional array, set to the number of times chains started in the given states
// res_N (output)
//     n x n matrix, filled with the sums of the number of transitions between states (i != j).  Diagonal entries are the exit-to-absorption counts
// workD (output)
//     3n^2 + 9n element workspace
// workI (output)
//     2n^2 element workspace
void LJMA_MHsample_Aslett(double *y, int *censored, int *m, double *pi, double *S, double *s, double *P, double *Pfull, double *Q, double *evals, double *Qinv_1, double *Qinv, int *n, int *iter, int *reverse, double *res_z, int *res_B, int *res_N, double *workD, int *workI) {
	LJMA_GetRNGstate();
	
	// Initialise output vars
	int i,j,k,l;
	for(i=0; i<*n; i++) {
		res_B[i] = 0;
		res_z[i] = 0.0;
		for(j=0; j<*n; j++) {
			res_N[i + j * *n] = 0;
		}
	}
	
	// Reverse storage needed?
	double *piR;
	if(*reverse == 0) {
		piR = NULL;
	} else {
		// Storage for stating state in reverse
		piR = workD; workD += *n;
	}
	
	// Create memory for current/proposal steps for MH
	double *c_z, *p_z, *tmpd;
	int c_B, p_B, *c_N, *p_N, c_pre, p_pre, *tmpi, tmpii;
	c_z = workD; workD += *n;
	p_z = workD; workD += *n;
	c_N = workI; workI += *n * *n;
	p_N = workI; workI += *n * *n;
	
	// Off we go
	double *y_p, U;
	int *censored_p;
	y_p = y;
	censored_p = censored;
	for(i=0; i<*m; i++) {
		//if(*reverse != 0) {
		//	LJMA_reversestart(y_p, pi, s, Q, evals, Qinv, n, piR, workD);
		//}
		// Get initial current step
		LJMA_samplechain(y_p, censored_p, pi, S, Q, evals, Qinv_1, P, Pfull, n, reverse, piR, c_z, &c_B, c_N, &c_pre, workD);
		//int cnt=0;
		while( (*reverse == 0 && s[c_pre] == 0.0) || (*reverse != 0 && pi[c_B] == 0.0) ) { // ensure the current is at least valid ... can't have a chain where the absorption was impossible from that state, or where we're starting from an impossible state
			LJMA_samplechain(y_p, censored_p, pi, S, Q, evals, Qinv_1, P, Pfull, n, reverse, piR, c_z, &c_B, c_N, &c_pre, workD); //cnt++;
		}
		//Rprintf("%d.", cnt);
		
		if( *censored_p == 0 ) { // only do MH step if not censored
			for(j=0; j<*iter; j++) {
				// i) get proposal
				LJMA_samplechain(y_p, censored_p, pi, S, Q, evals, Qinv_1, P, Pfull, n, reverse, piR, p_z, &p_B, p_N, &p_pre, workD);
				while( (*reverse == 0 && s[p_pre] == 0.0) || (*reverse != 0 && pi[p_B] == 0.0) ) {
					LJMA_samplechain(y_p, censored_p, pi, S, Q, evals, Qinv_1, P, Pfull, n, reverse, piR, p_z, &p_B, p_N, &p_pre, workD);
				}
				
				// ii) draw uniform
				U = runif(0.0, 1.0);
				
				// iii) select using accept ratio
				if( (*reverse == 0 && U < s[p_pre]/s[c_pre]) || (*reverse != 0 && U < pi[p_B]/pi[c_B]) ) {
					// for speed, if we move to proposal, just swap the memory that c and p point to!
					tmpd = c_z;
					c_z = p_z;
					p_z = tmpd;
					
					tmpii = c_B;
					c_B = p_B;
					p_B = tmpii;
					
					tmpi = c_N;
					c_N = p_N;
					p_N = tmpi;
					
					tmpii = c_pre;
					c_pre = p_pre;
					p_pre = tmpii;
				}
			}
		}
		
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
	
	LJMA_PutRNGstate();
}
