#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

#include "Simulate_AbsCTMC_gt_Hobolth_DCS.h"
#include "utility.h"

typedef struct HobPars_ {
	int lastj;
	int j;
	double prob;
	double *Q;
	double *evals;
	double *Qinv_b;
	double *S;
	double Pab;
	double *J;
	double y;
	double t;
	double u;
	int n;
} HobPars;
double HobCDF(double x, void *info) {
	HobPars *pars = (HobPars*) info;
	for(int i=0; i<pars->n; i++) {
		//if((pars->evals)[i]-(pars->S)[pars->lastj + pars->lastj * pars->n] == 0.0) {
		if(fabs(((pars->evals)[i]-(pars->S)[pars->lastj + pars->lastj * pars->n])/(pars->S)[pars->lastj + pars->lastj * pars->n]) < 1e-13) {
			(pars->J)[i] = x * exp((pars->evals)[i] * (pars->y-pars->t));
		} else {
			(pars->J)[i] = (exp((pars->evals)[i] * (pars->y-pars->t)) - exp((pars->y-pars->t-x) * (pars->evals)[i] + (pars->S)[pars->lastj + pars->lastj * pars->n] * x)) / ((pars->evals)[i]-(pars->S)[pars->lastj + pars->lastj * pars->n]);
		}
	}
	
	double tmp = 0.0;
	for(int i=0; i<pars->n; i++) {
		tmp += (pars->Q)[pars->j + i * pars->n] * (pars->J)[i] * (pars->Qinv_b)[i];
	}
	
	return(1 / pars->prob * (pars->S)[pars->lastj + pars->j * pars->n] / pars->Pab * tmp - pars->u);
}


//// Sample a chain conditional on absorbing after an observation, returning only sufficient statistics of the chain
// y (input)
//     absorption time to be conditioned on
// pi (input)
//     1 x n vector of initial state probabilities (*must* sum to 1)
// S (input)
//     n x n matrix of transition rates between non-absorbing states
// Q (input)
//     precomputed n x n matrix of eigenvectors of the original S matrix
// evals (input)
//     array size n of the n eigenvalues of the original S matrix
// Qinv_b (input)
//     precomputed n x 1 matrix Qinv %*% b, where b is vector with 1 for each exitable state
// b (input)
//     1 x n vector of 0/1s where 1 indicates an exit to be allowed
// s (input)
//     1 x n vector of exit rates
// n (input)
//     dimension
// res_z (output)
//     n dimensional array, filled with times spent in each state
// res_B (output)
//     integer, set to the state the sample chain started in
// res_N (output)
//     n x n matrix, filled with number of transitions between states (i != j).  The only non-zero diagonal was the exit-to-absorption state
// res_pre (output)
//     integer, set to the state chain was in immediately prior to absorption
// workD (output)
//     ?? element workspace
//
// return: sufficient statistics of the sample (z, B, N)
void LJMA_samplechain_Hobolth(double *y, double *pi, double *S, double *Q, double *evals, double *Qinv_b, double *b, double *s, int *n, double *res_z, int *res_B, int *res_N, int *res_pre, double *workD) {
	LJMA_GetRNGstate();
	
	int LOUIS = 0;
	// Initialise output vars
	*res_B = 0;
	for(int i=0; i<*n; i++) {
		res_z[i] = 0;
		for(int j=0; j<*n; j++) {
			res_N[i + j * *n] = 0;
		}
	}
	
	// Choose starting state ...
	double sofar = 0.0, target;
	// ... from pi
	target = runif(0.0, 1.0);
	*res_B = 0;
	while(sofar < target) {
		sofar += pi[(*res_B)++];
	}
	(*res_B)--;
	//Rprintf("%d - Start at: %d\n", *reverse, *res_B);
	

	HobPars pars;
	double (*HobCDFp)(double, void*);
	HobCDFp = &HobCDF;

	// Run through the chain
	int Maxit = 0;
	double Pab = 0.0, tmp = 0.0, p_sum, Tol = 0.0;
	double t = 0.0, jtime = 0.0;
	int j, lastj, i, k;
	double *J, *p;
	J = workD; workD += *n;
	p = workD; workD += *n;
	j = *res_B;
	while(t < *y) {
	//while(TRUE) {
		//lastt = t;
		lastj = j;
		
		// Commonly reused
		Pab = 0.0;
		for(i=0; i<*n; i++) {
			Pab += Q[j + i * *n] * exp(evals[i] * (*y-t)) * Qinv_b[i];
		}
		
		// First, we check if we are already in a state we can exit from ... if so, we might stay still
		if(b[j] > 0.0) {
//Rprintf("Could exit from %d, checking ...\n",j+1);
			if( runif(0.0, 1.0) < exp(S[j + j * *n]*(*y-t))/Pab ) {
				res_z[j] += (*y-t);
				res_N[j + j * *n] = 1;
				*res_pre = j; LOUIS = 1;
				break;
			}
		}
//Rprintf("Don't exit ...\n");
		
		// Ok, we're definitely moving somewhere then
		// Find J
		for(i=0; i<*n; i++) {
			//if(evals[i]-S[j + j * *n] == 0.0) {
			if( fabs((evals[i]-S[j + j * *n])/S[j + j * *n]) < 1e-13 ) {
				J[i] = (*y-t) * exp(evals[i] * (*y-t));
			} else {
				J[i] = (exp(evals[i] * (*y-t)) - exp(S[j + j * *n] * (*y-t))) / (evals[i]-S[j + j * *n]);
			}
		}
//for(i=0; i<*n; i++) Rprintf("J[%d] = %e, ", i+1, J[i]);
		
		// Find p_i
		p_sum = 0.0;
		for(i=0; i<*n; i++) {
			if(i == j) continue;
			
			tmp = 0.0;
			for(k=0; k<*n; k++) {
				tmp += Q[i + k * *n] * J[k] * Qinv_b[k];
			}
			
			p_sum += p[i] = S[j + i * *n] / Pab * tmp;
		}
		p[j] = 0.0;
//for(i=0; i<*n; i++) Rprintf("-> %d = %e, ", i+1, p[i]);
		
		// Make the jump
		sofar = 0.0;
		target = runif(0.0, p_sum);
		j = 0;
		while(sofar < target) {
			sofar += p[j++];
		}
		j--;
//Rprintf("\nChoose %d -> %d\n", lastj+1, j+1);
		
		// Get time advance
		pars.lastj = lastj;
		pars.j = j;
		pars.prob = p[j];
		pars.Q = Q;
		pars.evals = evals;
		pars.Qinv_b = Qinv_b;
		pars.S = S;
		pars.Pab = Pab;
		pars.J = J;
		pars.y = *y;
		pars.t = t;
		pars.u = runif(0.0, 1.0);
		pars.n = *n;
		
		Tol = 0.0;
		Maxit = 1000;
		jtime = R_zeroin2(0.0, *y-t, -pars.u, 1.0-pars.u, HobCDFp, &pars, &Tol, &Maxit);
		if(Maxit == -1) {
			Rprintf("\nWARNING: CDF root finder didn't converge\ny-t=%e, u=%e, Tol=%e, jtime=%e, prob=%e\n", *y-t, pars.u, Tol, jtime, Pab);
			//Rprintf("\n\nERROR ERROR, ABORT ABORT: t=%.20e\ny=%.20e\ny-t+jtime=%.20e\ny-t=%.20e\nu=%.20e\njtime=%.20e\nTol=%.20e\nIt=%d\nL_bound=%.20e\nR_bound=%.20e\n\n", pars.t, pars.y, t, pars.y-pars.t, pars.u, jtime, Tol, Maxit, HobCDF(0, &pars), HobCDF(pars.y-pars.t, &pars));
			for(int i=0; i<*n; i++) {
				Rprintf("*%e*", ((pars.evals)[i]-(pars.S)[pars.lastj + pars.lastj * pars.n]));
			}
			for(double ii=0.0; ii<=1.0; ii+=0.1) {
				Rprintf("x <- c(x, %.15e)\ny <- c(y, %.15e)\n", (pars.y-pars.t)*ii, HobCDF((pars.y-pars.t)*ii, &pars));
			}
			Rprintf("\n\n");
		}
//Rprintf("Jump time %e\n", jtime);
		
		// Sometimes the time left is so small that under finite precision almost any time we add will cause y-t=0, which should never happen here ... only exit by check above
		while(t+jtime >= *y) {
			jtime = jtime/2;
		}
		
		// Update sufficient statistics
		res_N[lastj + j * *n]++;
		res_z[lastj] += jtime;
		
		// Advance time, continue
		t += jtime;
	}
	
	if(LOUIS==0) {
		Rprintf("\n\nERROR ERROR, ABORT ABORT: t=%.20e\ny=%.20e\ny-t+jtime=%.20e\ny-t=%.20e\nu=%.20e\njtime=%.20e\nTol=%.20e\nIt=%d\nL_bound=%.20e\nR_bound=%.20e\n\n", pars.t, pars.y, t, pars.y-pars.t, pars.u, jtime, Tol, Maxit, HobCDF(0, &pars), HobCDF(pars.y-pars.t, &pars));
		for(double ii=0.0; ii<=1.0; ii+=0.1) {
			Rprintf("x <- c(x, %.15e)\ny <- c(y, %.15e)\n", (pars.y-pars.t)*ii, HobCDF((pars.y-pars.t)*ii, &pars));
		}
		Rprintf("\n\n");
	}
	
	LJMA_GUI();
	LJMA_PutRNGstate();
}




//// Metropolis-Hastings sample a chain using directly conditional sampling, returning only sufficient statistics of the chain
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
void LJMA_MHsample_Hobolth(double *y, int *censored, int *m, double *pi, double *S, double *s, double *Q, double *evals, double *Qinv_b, double *b, double *Qinv, int *n, int *iter, double *res_z, int *res_B, int *res_N, double *workD, int *workI) {
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
		// Get initial current step
		if(*censored_p) Rprintf("\nWARNING: Hobolth does not support censoring\n");
		LJMA_samplechain_Hobolth(y_p, pi, S, Q, evals, Qinv_b, b, s, n, c_z, &c_B, c_N, &c_pre, workD);
		//int cnt=0;
		while( s[c_pre] == 0.0 ) { // ensure the current is at least valid ... can't have a chain where the absorption was impossible from that state, or where we're starting from an impossible state
			LJMA_samplechain_Hobolth(y_p, pi, S, Q, evals, Qinv_b, b, s, n, c_z, &c_B, c_N, &c_pre, workD); //cnt++;
		}
		//Rprintf("%d.", cnt);
		
		//Rprintf("y=%e\nz=%e %e %e\nN=\n", *y_p, c_z[0], c_z[1], c_z[2]);
		//for(int i=0; i<*n; i++) { for(int j=0; j<*n; j++) { Rprintf("%2d ", c_N[i + j**n]); } Rprintf("\n"); }
		//Rprintf("\n\n");

		if( *censored_p == 0 ) { // only do MH step if not censored
			for(j=0; j<*iter; j++) {
				// i) get proposal
				LJMA_samplechain_Hobolth(y_p, pi, S, Q, evals, Qinv_b, b, s, n, p_z, &p_B, p_N, &p_pre, workD);
				while( s[c_pre] == 0.0 ) { // ensure the current is at least valid ... can't have a chain where the absorption was impossible from that state, or where we're starting from an impossible state
					LJMA_samplechain_Hobolth(y_p, pi, S, Q, evals, Qinv_b, b, s, n, p_z, &p_B, p_N, &p_pre, workD); //cnt++;
				}
				
				// ii) draw uniform
				U = runif(0.0, 1.0);
				
				// iii) select using accept ratio
				if( U < s[p_pre]/s[c_pre] ) {
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

