#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "Simulate_AbsCTMC_gt_Bladt_MHRS.h"
#include "Simulate_AbsCTMC_eq_Bladt_MHRS.h"

//// Metropolis-Hastings + Rejection Sampling sample a chain conditional on absorbing after an observation, returning only sufficient statistics of the chain
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
// Pfull (input)
//     n x (n+1) matrix of embedded transition probabilities
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
//     ??? element workspace
// workI (output)
//     ??? element workspace
//
void LJMA_MHsample_Bladt(double *y, int *censored, int *m, double *pi, double *S, double *s, double *Pfull, int *n, int *iter, double *res_z, int *res_B, int *res_N, double *workD, int *workI) {
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
	
	// Create memory for current/proposal steps for MH
	double *c_z, *p_z, *tmpd, absorbed;
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
		// Get initial current step     void LJMA_samplechain_Bladt(double *y, int *censored, double *pi, double *S, int *n, double *absorbed, double *res_z, int *res_B, int *res_N, int *res_pre)
		LJMA_samplechain_Bladt(y_p, censored_p, pi, S, Pfull, n, &absorbed, c_z, &c_B, c_N, &c_pre, workD, workI);
		while(s[c_pre] == 0) { // ensure the current is at least valid ... can't have a chain where the absorption was impossible from that state
			LJMA_samplechain_Bladt(y_p, censored_p, pi, S, Pfull, n, &absorbed, c_z, &c_B, c_N, &c_pre, workD, workI);
		}
		
		if( *censored_p == 0 ) { // We don't need to do the MH-accept step if this was a censored sample
			for(j=0; j<*iter; j++) {
				// i) get proposal
				LJMA_samplechain_Bladt(y_p, censored_p, pi, S, Pfull, n, &absorbed, p_z, &p_B, p_N, &p_pre, workD, workI);
				while(s[p_pre] == 0) {
					LJMA_samplechain_Bladt(y_p, censored_p, pi, S, Pfull, n, &absorbed, p_z, &p_B, p_N, &p_pre, workD, workI);
				}
				
				// ii) draw uniform
				U = runif(0.0, 1.0);
				
				// iii) select using accept ratio
				if(U < s[p_pre]/s[c_pre]) {
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
