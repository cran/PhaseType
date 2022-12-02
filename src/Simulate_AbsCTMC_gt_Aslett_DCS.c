#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "arms.h"
#include "utility.h"
#include "Simulate_AbsCTMC_gt_Aslett_DCS.h"

//// Optimally compute 1 - pi %*% Q %*% exp(x * evals) %*% Qinv %*% 1
//// In other words, compute upper tail PHT cdf value for optimal case where pi and S never change, only x
// x (input)
//     scalar value
// pi_Q (input)
//     precomputed 1 x n matrix pi %*% Q of the original S matrix
// evals (input)
//     array size n of the n eigenvalues of the original S matrix
// Qinv_1 (input)
//     precomputed n x 1 matrix Qinv %*% 1
// n (input)
//     dimension
// res (output)
//     stores result here too
//
// return: upper tail PHT cdf value
void LJMA_phtcdf_opt(double *x, double *pi_Q, double *evals, double *Qinv_1, int *n, double *res) {
	if(*x > 0) {
		double *pi_Q_p, *evals_p, *Qinv_1_p;
		pi_Q_p = pi_Q;
		evals_p = evals;
		Qinv_1_p = Qinv_1;

		int i;
		double result = 0.0;
		for(i=0; i<*n; i++) {
			result += *(pi_Q_p++) * exp(*x * *(evals_p++)) * *(Qinv_1_p++);
		}

		*res = result;
	} else {
		*res = 1.0;
	}
}


//// Compute 1 - pi %*% Q %*% exp(x * evals) %*% Qinv %*% 1
//// In other words, compute upper tail PHT cdf value for S never changes, only x
// x (input)
//     scalar value of time to evaluate CDF at
// pi (input)
//     1 x n vector of initial state probabilities
// Q (input)
//     precomputed n x n matrix of eigenvectors of the original S matrix
// evals (input)
//     array size n of the n eigenvalues of the original S matrix
// Qinv_1 (input)
//     precomputed n x 1 matrix Qinv %*% 1
// n (input)
//     dimension
// res (output)
//     stores result here too
// workD (output)
//     n element workspace
//
// return: upper tail PHT cdf value
void LJMA_phtcdf(double *x, double *pi, double *Q, double *evals, double *Qinv_1, int *n, double *res, double *workD) {
	if(*x > 0) {
		// First, compute pi %*% Q, then we can do the usual
		double *pi_Q;
		pi_Q = workD;

		char trans = 'T';
		double alpha = 1.0, beta = 0.0;
		int incx = 1, incy = 1;
		F77_CALL(dgemv)(&trans, n, n, &alpha, Q, n, pi, &incx, &beta, pi_Q, &incy FCONE);

		LJMA_phtcdf_opt(x, pi_Q, evals, Qinv_1, n, res);
	} else {
		*res = 1.0;
	}
}


//// Compute conditional jump density given absorption after a certain time
// d (input)
//     scalar value at which to evaluate density
// tnow (input)
//     the current time at which the conditional jump is being computed
// jnow (input)
//     the current state at which the conditional jump is being computed
// y (input)
//     absorption time to be conditioned on
// S (input)
//     n x n matrix of transition rates between non-absorbing states
// Q (input)
//     precomputed n x n matrix of eigenvectors of the original S matrix
// evals (input)
//     array size n of the n eigenvalues of the original S matrix
// Qinv_1 (input)
//     precomputed n x 1 matrix Qinv %*% 1
// P (input)
//     n x n matrix of embedded chain transition probabilities (excl. absorbing move)
// n (input)
//     dimension
// res (output)
//     stores result here too
// workD (output)
//     3n element workspace
//
// return: conditional jump density
void LJMA_condjumpdens(double *d, double *tnow, int *jnow, double *y, double *S, double *Q, double *evals, double *Qinv_1, double *P, int *n, double *res, double *workD) {
	// First, set up the pi values for the numerator and denominator PHT dist values
	double *pi1, *pi2;
	pi1 = workD; workD += *n;
	pi2 = workD; workD += *n;

	for(int i=0; i<*n; i++) {
		pi1[i] = P[*jnow + i * *n];
		pi2[i] = 0.0;
	}
	pi2[*jnow] = 1.0;

	// Now compute the conditional jump density
	double r1, x1; //, x2;
	x1 = *y - *tnow - *d;
	//x2 = *y - *tnow;
	////**** Could wring more performance out of this ... r2 doesn't change over all rejection samples for a given jump
	if(x1 > 0) LJMA_phtcdf(&x1, pi1, Q, evals, Qinv_1, n, &r1, workD); else r1 = 1;
	//if(x2 > 0) LJMA_phtcdf(&x2, pi2, Q, evals, Qinv_1, n, &r2, workD); else r2 = 1; // DON'T EVEN NEED TO CALC THIS!! JUST NEED UP TO A CONST.
	//*res = (r1 * dexp(*d, -1.0/S[*jnow + *jnow * *n], 0)) / r2;
	*res = log(r1) + dexp(*d, -1.0/S[*jnow + *jnow * *n], 1); // - log(r2); // LEAVE DENOMINATOR OUT OF IT!! JUST NEED UP TO A CONST.
}


//// Supporting function for adaptive rejection sampling method
typedef struct condjumpdens_pars_ {
	double *tnow;
	int *jnow;
	double *y;
	double *S;
	double *Q;
	double *evals;
	double *Qinv_1;
	double *P;
	int *n;
	double *workD;
} condjumpdens_pars;
double LJMA_condjumpdens_ars(double d, void *par) {
	double res;
	condjumpdens_pars *pars;
	pars = (condjumpdens_pars *) par;
	LJMA_condjumpdens(&d, pars->tnow, pars->jnow, pars->y, pars->S, pars->Q, pars->evals, pars->Qinv_1, pars->P, pars->n, &res, pars->workD);
	/*Rprintf("(%lf, %lf)\n", d, log(res));*/
	//return(log(res));
	return(res);
}


//// Sample from conditional jump density given absorption after a certain time, by adaptive rejection sampling
// tnow (input)
//     the current time at which the conditional jump is being computed
// jnow (input)
//     the current state at which the conditional jump is being computed
// y (input)
//     absorption time to be conditioned on
// S (input)
//     n x n matrix of transition rates between non-absorbing states
// Q (input)
//     precomputed n x n matrix of eigenvectors of the original S matrix
// evals (input)
//     array size n of the n eigenvalues of the original S matrix
// Qinv_1 (input)
//     precomputed n x 1 matrix Qinv %*% 1
// P (input)
//     n x n matrix of embedded chain transition probabilities (excl. absorbing move)
/// n (input)
//     dimension
// res (output)
//     stores result here too
// workD (output)
//     3n element workspace
//
// return: a sample from the conditional jump density
void LJMA_condjump_r_ars(double *tnow, int *jnow, double *y, double *S, double *Q, double *evals, double *Qinv_1, double *P, int *n, double *res, double *workD) {
	LJMA_GetRNGstate();

	if(*tnow >= *y) {
		*res = rexp(1.0/-S[*jnow+*jnow*(*n)]);
		LJMA_PutRNGstate();
		return;
	}

	double x = *y-*tnow, *pi = workD, denom;
	for(int i=0; i<*n; i++) {
	  pi[i] = 0.0;
	}
	pi[*jnow] = 1.0;
	LJMA_phtcdf(&x, pi, Q, evals, Qinv_1, n, &denom, workD + *n);
	//	Rprintf("%lf, ", exp(S[*jnow+*jnow*(*n)] * (*y - *tnow))/denom);
	if(*tnow < *y && runif(0.0,1.0) < exp(S[*jnow+*jnow*(*n)] * (*y - *tnow))/denom) {
		*res = *y - *tnow + rexp(1.0/-S[*jnow+*jnow*(*n)]);
		LJMA_PutRNGstate();
		return;
	}
	if(*tnow > *y) Rprintf("ALERT!\n");

	condjumpdens_pars pars;
	pars.tnow = tnow;
	pars.jnow = jnow;
	pars.y = y;
	pars.S = S;
	pars.Q = Q;
	pars.evals = evals;
	pars.Qinv_1 = Qinv_1;
	pars.P = P;
	pars.n = n;
	pars.workD = workD;

	double xinit[4], xl, xr, convex, xprev, xsamp, qcent, xcent;
	int info, ninit, npoint, dometrop, nsamp, ncent, neval;

	//xinit[3] = *y;
	//while(!isfinite(LJMA_condjumpdens_ars(xinit[3], &pars))) {
	//	xinit[3] = log(xinit[3]);
	//}
	//xinit[0] = xinit[3]/4.0;
	xinit[0] = (*y - *tnow)/1e6;
	xinit[1] = (*y - *tnow)/3.0;
	xinit[2] = xinit[1]*2.0;
	xinit[3] = *y - *tnow - xinit[0];
	ninit = 4;

	xl = 0.0;
	xr = *y - *tnow;

	convex = 1.0;
	npoint = 100;
	dometrop = 1;
	xprev = 0.0;
	xsamp = 0.0;
	nsamp = 1;
	qcent = 0.0;
	xcent = 0.0;
	ncent = 0;
	neval = 0;

	double (*myfunc)(double, void*);
	myfunc = &LJMA_condjumpdens_ars;

	info = arms(xinit, ninit, &xl, &xr, myfunc, &pars, &convex, npoint, dometrop, &xprev, &xsamp, nsamp, &qcent, &xcent, ncent, &neval);
	if(info != 0) {
		Rprintf("Error (LJMA_condjump_r_ars 01): ARS failed, code=%d\n", info);
	}
	//Rprintf("%lf - ", xsamp);
	if(xsamp > *y - *tnow) Rprintf("ALERT!!\n");

	*res = xsamp;

	LJMA_PutRNGstate();
}


//// Sample a chain conditional on absorbing after an observation, returning only sufficient statistics of the chain
// y (input)
//     absorption time to be conditioned on
// censored (input)
//    0=not censored, 1=censored
// pi (input)
//     1 x n vector of initial state probabilities (*must* sum to 1)
// S (input)
//     n x n matrix of transition rates between non-absorbing states
// Q (input)
//     precomputed n x n matrix of eigenvectors of the original S matrix
// evals (input)
//     array size n of the n eigenvalues of the original S matrix
// Qinv_1 (input)
//     precomputed n x 1 matrix Qinv %*% 1
// P (input)
//     n x n matrix of embedded chain transition probabilities (excl. absorbing move)
// Pfull (input)
//     n x (n+1) matrix of embedded chain transition probabilities, incluing option of absorption
// n (input)
//     dimension
// reverse (input)
//     whether to sample in reverse
// piR (input)
//     "initial" state distribution for reverse sampling (ignored unless reverse != 0)
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
//
void LJMA_samplechain(double *y, int *censored, double *pi, double *S, double *Q, double *evals, double *Qinv_1, double *P, double *Pfull, int *n, int *reverse, double *piR, double *res_z, int *res_B, int *res_N, int *res_pre, double *workD) {
	LJMA_GetRNGstate();

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
	if(*reverse == 0) {
		// ... from pi
		target = runif(0.0, 1.0);
		*res_B = 0;
		while(sofar < target) {
			sofar += pi[(*res_B)++];
		}
		(*res_B)--;
	} else {
		// ... from piR
		target = runif(0.0, 1.0);
		*res_pre = 0;
		while(sofar < target) {
			sofar += piR[(*res_pre)++];
		}
		(*res_pre)--;
	}
	//Rprintf("%d - Start at: %d\n", *reverse, *res_B);

	// Run through the chain
	double t = 0.0, lastt = 0.0;
	int j, lastj = 0, i;
	double d, x1, r1, r2, *pi1, *pi2;
	pi1 = workD; workD += *n;
	pi2 = workD; workD += *n;
	if(*reverse == 0) { j = *res_B; } else { j = *res_pre; }
	while(t < *y || *censored) {
		lastt = t;
		lastj = j;

		// Make the time jump
		/*LJMA_condjump_r(&t, &j, y, S, Q, evals, Qinv_1, n, &d);*/
		LJMA_condjump_r_ars(&t, &j, y, S, Q, evals, Qinv_1, P, n, &d, workD); LJMA_counter++;
		//if(*censored == 1) { Rprintf("- %d/%.2lf/%.2lf/%.2lf/", j, t, *y, d); }
		t += d;

		// Make the state jump
		target = runif(0.0, 1.0); //if(*censored == 1) { Rprintf("%.2lf(", target); }
		sofar = 0.0;
		j = 0;
		if(t < *y) {
			x1 = *y - t;
			for(i=0; i<*n; i++) {
				pi2[i] = P[lastj + i * *n];
			}
			LJMA_phtcdf(&x1, pi2, Q, evals, Qinv_1, n, &r2, workD);
			while(sofar < target) {
				if(P[lastj + j * *n] == 0.0) { j++; continue; } // Don't bother calcing if impossible move!
				for(i=0; i<*n; i++) {
					pi1[i] = 0.0;
				}
				pi1[j] = 1.0;
				LJMA_phtcdf(&x1, pi1, Q, evals, Qinv_1, n, &r1, workD);
				sofar += r1*P[lastj + (j++) * *n]/r2;
				//Rprintf("%d -> %d: %lf\n", lastj, j-1, r1*P[lastj + (j-1) * *n]/r2);
				// OLD AND WRONG: sofar += P[lastj + j++ * *n];
			}
		} else {
			while(sofar < target) {
				sofar += Pfull[lastj + j++ * *n];
				//Rprintf("%d*%.2lf,", j-1, sofar);
			}
		}
		//Rprintf(") -");
		j--;

		if(j == *n) {
			break;
		} else if(*reverse == 0 && (t < *y || *censored)) {
			res_z[lastj] += t - lastt;
			res_N[lastj + j * *n]++;
		} else if(*reverse != 0 && (t < *y || *censored)) {
			res_z[lastj] += t - lastt;
			res_N[j + lastj * *n]++;
		}
	}

	if(*censored == 0) { if(*reverse == 0) {	res_z[lastj] += *y - lastt; *res_pre = lastj; } else { res_z[lastj] += *y - lastt; *res_B = lastj; } }
	else { if(*reverse == 0) {	res_z[lastj] += t - lastt; *res_pre = lastj; } else { res_z[lastj] += t - lastt; *res_B = lastj; } }
	res_N[*res_pre + *res_pre * *n]++; // records on the diagonal what state we jumped to absorption from

	/*if(*censored == 1) {
	 Rprintf("\n\nB = %d\n", *res_B);
	 Rprintf("z = ");
	 for(int i=0; i<*n; i++) {
	 Rprintf("%lf ", res_z[i]);
	 }
	 Rprintf("\nN = \n");
	 for(int i=0; i<*n; i++) {
	 for(int j=0; j<*n; j++) {
	 Rprintf("%d ", res_N[i + j * *n]);
	 }
	 Rprintf("\n");
	 }
	 Rprintf("\nPfull = \n");
	 for(int i=0; i<*n; i++) {
	 for(int j=0; j<*n+1; j++) {
	 Rprintf("%.2lf ", Pfull[i + j * *n]);
	 }
	 Rprintf("\n");
	 }
	 }*/

	LJMA_GUI();
	LJMA_PutRNGstate();
}


//// Compute the "starting" state distribution for reverse sampling
// t (input)
//     time at which to compute "starting" state
// pi (input)
//     1 x n vector of initial state probabilities (*must* sum to 1)
// S (input)
//     n x n matrix of transition rates between non-absorbing states
// s (input)
//     1 x n vector of exit rates
// n (input)
//     dimension
// res_pi (output)
//     1 x n vector of "initial" state probabilities for reverse sampling (normalised to sum to 1)
// workD (output)
//     n^2+n element workspace
void LJMA_reversestart(double *t, double *pi, double *s, double *Q, double *evals, double *Qinv, int *n, double *res_pi, double *workD) {
	// R code for computing reverse starting distribution
	//piR2 <- vector("numeric", n)
	//piR2norm <- pi %*% e$vectors %*% diag(exp(e$values*t),n,n) %*% solve(e$vectors) %*% s
	//for(i in 1:n) {
	//	piR2[i] <- pi %*% e$vectors %*% diag(exp(e$values*t),n,n) %*% solve(e$vectors) %*% {tmp<-rep(0,n); tmp[i]<-s[i]; tmp} / piR2norm
	//}
	//print(piR2)
	double *a, denom = 0.0, *b;
	a = workD; workD += *n; b = workD; workD += *n;
	char trans = 'T'; double alpha = 1.0, beta = 0.0; int inc = 1;
	F77_CALL(dgemv)(&trans, n, n, &alpha, Q, n, pi, &inc, &beta, a, &inc FCONE); // pi %*% e$vectors -> a
	// a %*% diag(exp(e$values*t),n,n) -> a
	for(int i=0; i<*n; i++) {
		a[i] = a[i] * exp(evals[i] * *t);
	}
	F77_CALL(dgemv)(&trans, n, n, &alpha, Qinv, n, a, &inc, &beta, b, &inc FCONE); // a %*% solve(e$vectors)
	// then find denominator
	for(int i=0; i<*n; i++) {
		denom += b[i] * s[i];
	}
	for(int i=0; i<*n; i++) {
		res_pi[i] = b[i] * s[i] / denom;
	}
}


//// Compute the stationary distribution of a CTMC
// S (input)
//     n x n matrix of transition rates - NOTE: diagonals are ignored and recomputed from off-diagonals, otherwise incorrect dist calculated
// n (input)
//     dimension
// res_pi (output)
//     1 x n vector of "initial" state probabilities for reverse sampling (normalised to sum to 1)
// workD (output)
//     n^2 + n element workspace
// workI (output)
//     n element workspace
void LJMA_stationary(double *S, int *n, double *res_pi, double *workD, int *workI) {
	double *Q;
	Q = workD; workD += *n * *n;
	memcpy(Q, S, *n * *n * sizeof(*S));

	// Q + E
	for(int i=0; i < *n; i++) {
		Q[i + i * *n] = 0.0;
		for(int j=0; j < *n; j++) {
			if(i!=j) {
				Q[i + i * *n] -= Q[i + j * *n]++;
			}
		}
		Q[i + i * *n]++;
	}

	/*Rprintf("\n\nQ+E: ");
	 for(int i=0; i<*n; i++) {
	 for(int j=0; j<*n; j++) {
	 Rprintf("%lf ", Q[i + j * *n]);
	 }
	 Rprintf("\n");
	 }
	 Rprintf("\n\n");*/


	// (Q+E)^{-1}
	LJMA_inverse(Q, n, workI);

	// e.(Q+E)^{-1}
	double *e;
	e = workD; workD += *n;
	for(int i=0; i<*n; i++) { e[i] = 1.0; }
	char trans = 'T'; double alpha = 1.0, beta = 0.0; int inc = 1;
	F77_CALL(dgemv)(&trans, n, n, &alpha, Q, n, e, &inc, &beta, res_pi, &inc FCONE);
}
