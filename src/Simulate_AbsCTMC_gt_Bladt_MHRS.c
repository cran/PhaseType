#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "Simulate_AbsCTMC_gt_Bladt_MHRS.h"

//// Rejection Sampling sample a chain conditional on absorbing after an observation, returning only sufficient statistics of the chain
// y (input)
//     absorption time to be reached by rejection sampling
// censored (input)
//     1/0 for is/not a censored observation
// pi (input)
//     1 x n vector of initial state probabilities (*must* sum to 1)
// S (input)
//     n x n matrix of transition rates between non-absorbing states
// Pfull (input)
//     n x (n+1) matrix of embedded chain transition probabilities, incluing option of absorption
// n (input)
//     dimension
// absorbed (output)
//     time at which absorption took place
// res_z (output)
//     n dimensional array, filled with times spent in each state
// res_B (output)
//     integer, set to the state the sample chain started in
// res_N (output)
//     n x n matrix, filled with number of transitions between states (i != j).  The only non-zero diagonal was the exit-to-absorption state
// res_pre (output)
//     integer, set to the state chain was in immediately prior to absorption
// workD (output)
//     ??? element workspace
// workI (output)
//     ??? element workspace
//
void LJMA_samplechain_Bladt(double *y, int *censored, double *pi, double *S, double *Pfull, int *n, double *absorbed, double *res_z, int *res_B, int *res_N, int *res_pre, double *workD, int *workI) {
	LJMA_GetRNGstate();
	
	double t = 0.0; // Track time
	
	int B2 = 0, *N2, pre2; // These store the return values for proposed jumps ...
	double *z2; // ... since rejection sampling means we don't know when we'll get the run we return
	
	int lastj = 0, j;
	double lastt = 0.0, *p, sofar=0.0, target=0.0;
	p = workD; workD += *n + 1;
	
	N2 = workI; workI += *n * *n;
	z2 = workD; workD += *n;
	
	while(t < *y) {
		t = 0;
		
		// Initialise output vars
		for(int i=0; i<*n * *n; i++) N2[i] = 0;
		for(int i=0; i<*n; i++) z2[i] = 0.0;
		
		// Choose starting state ...
		target = runif(0.0, 1.0);
		sofar = 0.0;
		// ... from pi
		B2 = 0;
		while(sofar < target && B2 <= *n) {
			sofar += pi[B2++];
		}
		B2--;
		
#ifdef DEBUG
		Rprintf("====\nGenerated starting step: %d\n", B2);
#endif
		
		// Proceed through to absorption
		j = B2;
		lastt = t;
		lastj = j;
		//		while(j < *n) {
		while((t < *y && j < *n) || (*censored && j < *n)) {
			// Update time with time of next jump
#ifdef DEBUG
			Rprintf("Generating new jump time, exponential mean %lf\n", 1.0/-S[j+j*(*n)]);
#endif
			t = t + rexp(1.0/-S[j + j*(*n)]); LJMA_counter++;
#ifdef DEBUG
			Rprintf("New jump after %lf @ %lf\n", t-lastt, t);
#endif
			
			// Allocate jump probabilities
			for(int i=0; i < *n+1; i++) {
				p[i] = Pfull[j+i*(*n)];
			}
#ifdef DEBUG
			Rprintf("Jump probabilities =");
			for(int i=0; i < *n+1; i++) {
				Rprintf(" %lf ", p[i]);
			}
			Rprintf("\n");
#endif
			
			// Generate jump
			target = runif(0.0, 1.0);
			j = 0;
			sofar = 0.0;
			while(sofar < target && j <= *n + 1) {
				sofar += p[j++];
			}
			j--;
			
#ifdef DEBUG
			Rprintf("Jump chosen = %d\n", j);
#endif
			
			// Record what we need, if we're not past x yet (depending on whether censored or not)
			if((t < *y && j < *n) || (*censored && j < *n)) {
				z2[lastj] += t - lastt;
				N2[lastj + j * *n]++;
				if(j < *n) {
					lastj = j;
					lastt = t;
				}
			}
		}
		LJMA_GUI();
	}
	
#ifdef DEBUG
	Rprintf("Exited while, rejection sample found\n");
#endif
	
	// Out here means we have our chain ... add absorbing step
	/**res_pre = lastj;
	 N2[lastj+lastj*(*n)]++;
	 if(!*censored) {
	 z2[lastj] += (*y - lastt);
	 } else {
	 z2[lastj] += (t - lastt);
	 }*/
	if(*censored == 0) { z2[lastj] += *y - lastt; pre2 = lastj; }
	else { z2[lastj] += t - lastt; pre2 = lastj; }
	N2[pre2 + pre2 * *n]++; // records on the diagonal what state we jumped to absorption from
	
	
#ifdef DEBUG
	Rprintf("About to copy results for return ...\n");
#endif
	
	// Now copy our proposal to the return result as we are done and have made rejection sampling selection
	*absorbed = t;
	*res_B = B2;
	*res_pre = pre2;
	for(int i=0; i < *n; i++) {
		res_z[i] = z2[i];
		for(int j=0; j < *n; j++) {
			res_N[i+j*(*n)] = N2[i+j*(*n)];
		}
	}
	
#ifdef DEBUG
	Rprintf("Results copied, tidy and return ...\n");
#endif
	
	LJMA_PutRNGstate();
}
