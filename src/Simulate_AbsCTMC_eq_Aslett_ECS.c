#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include "utility.h"
#include "Simulate_AbsCTMC_eq_Aslett_ECS.h"
#include "Simulate_AbsCTMC_gt_Aslett_DCS.h"


/* R VERSION
 moveDens <- function(j, S, s, y, t, d) {
   dens <- S[j,]
   dens[j] <- 0
   dens <- dens * expm(S*(y-t-d)) %*% s
   
   dens <- dens/sum(dens)
   dens
 }
 */
void LJMA_moveMass(double *p, double y_t_d, int j, double *P, double *Q, double *evals, double *Qinv_s, int n, double *workD) {
	double *tmp;
	tmp = workD;
	for(int i=0; i<n; i++) {
		tmp[i] = exp(evals[i]*y_t_d) * Qinv_s[i];
	}
	
	char transN = 'N';
	double oneD = 1.0, zeroD = 0.0;
	int oneI = 1.0;
	F77_CALL(dgemv)(&transN, &n, &n, &oneD, Q, &n, tmp, &oneI, &zeroD, p, &oneI);
	
	double sum = 0.0;
	for(int i=0; i<n; i++) {
		sum += p[i] = p[i] * P[j + i * n];
	}
	for(int i=0; i<n; i++) {
		p[i] = p[i]/sum;
		//Rprintf("Prob -> %d = %lf\n", i+1, p[i]);
	}
}

// Dot product to be root found on
typedef struct dotProd_pars_ {
	double U;
	int j;
	double *PjQSjjLinv;
	double *evals_Sjj;
	double *diagConst;
	double *Ly_t;
	double *evals;
	double *Qinv_s;
	int n;
	double xi;
} dotProd_pars;
/* R VERSION
 dotProdSolv <- function(x,z,j,Q,L,s,P,y,t) {
   D1 <- (1/(rep(-S[j,j],3)+L))*exp(L*(y-t))
   Dd <- (1/(rep(-S[j,j],3)+L))*exp(L*(y-t)-x*(L+rep(-S[j,j],3)))
   D2 <- (1/(rep(-S[j,j],3)+L))*exp(-rep((y-t)*(-S[j,j]),3))
   
   V <- Q %*% diag((D1*z-D2*z),3,3) %*% solve(Q)
   Ud <- Q %*% diag((D1-Dd),3,3) %*% solve(Q)
   
   P[j,]%*%Q%*%diag(D1*(z-1)-D2*z+Dd,3)%*%solve(Q)%*%s # It should be this
   #(P[j,] %*% Ud %*% s) / (P[j,] %*% Q %*% diag(D1-D2,3,3) %*% solve(Q) %*% s) - z # full CDF way
 }
*/
/* ORIGINAL VERSION: */
double LJMA_dotProd(double d, void *par) {
	double res = 0.0;
	//dotProd_pars *pars;
	//pars = (dotProd_pars *) par;
	
	int offset = ((dotProd_pars *)par)->j * ((dotProd_pars *)par)->n;
	for(int i=0; i < ((dotProd_pars *)par)->n; i++) {
		if( ((dotProd_pars *)par)->evals_Sjj[offset] == 0.0 ) {
			res += ((dotProd_pars *)par)->PjQSjjLinv[offset] * ( ((dotProd_pars *)par)->diagConst[i] + d * exp(((dotProd_pars *)par)->Ly_t[i]) ) * ((dotProd_pars *)par)->Qinv_s[i];
		} else {
			res += ((dotProd_pars *)par)->PjQSjjLinv[offset] * ( ((dotProd_pars *)par)->diagConst[i] + exp(((dotProd_pars *)par)->Ly_t[i] -d * ((dotProd_pars *)par)->evals_Sjj[offset]) ) * ((dotProd_pars *)par)->Qinv_s[i];
		}
		//Rprintf("d=%lf, inner[%d] = %e, exp(inner) = %e, res = %e\n", d, i, ((dotProd_pars *)par)->Ly_t[i] -d * ((dotProd_pars *)par)->evals_Sjj[offset], exp(((dotProd_pars *)par)->Ly_t[i] -d * ((dotProd_pars *)par)->evals_Sjj[offset]), ((dotProd_pars *)par)->PjQSjjLinv[offset] * ( ((dotProd_pars *)par)->diagConst[i] + exp(((dotProd_pars *)par)->Ly_t[i] -d * ((dotProd_pars *)par)->evals_Sjj[offset]) ) * ((dotProd_pars *)par)->Qinv_s[i]);
		offset++;
	}
	
	return(res);
}
/*double LJMA_dotProd(double d, void *par) {
	double res = 0.0;
	dotProd_pars *pars;
	pars = (dotProd_pars *) par;
	
	int j = pars->j, n = pars->n;
	for(int i=0; i < pars->n; i++) {
		//  res += pars->PjQSjjLinv[j + i * n] * ( pars->diagConst[i] + pars->eLy_t[i] * exp( -d*(pars->evals[i] + pars->xi)) ) * pars->Qinv_s[i];
		res += pars->PjQSjjLinv[j + i * n] * ( pars->diagConst[i] + exp(pars->Ly_t[i] -d*pars->evals_Sjj[i + j * n]) ) * pars->Qinv_s[i];
	}
	
	return(res);
}*/
// Same function, callable from R
void LJMA_dotProd_R(double *d, double *U, int *j, double *PjQSjjLinv, double *diagConst, double *Ly_t, double *evals, double *Qinv_s, int *n, double *xi, double *res) {
	dotProd_pars pars;
	pars.U = *U;
	pars.j = *j;
	pars.PjQSjjLinv = PjQSjjLinv;
	pars.diagConst = diagConst;
	pars.Ly_t = Ly_t;
	pars.evals = evals;
	pars.Qinv_s = Qinv_s;
	pars.n = *n;
	pars.xi = *xi;
	
	*res = LJMA_dotProd(*d, &pars);
}



// {exp(S[j,j]*(y-t)) * s[j]}/{e_j %*% expm(S*(y-t)) %*% s}
double LJMA_probAbsorb(double y_t, int j, double *S, double *Q, double *evals, double *Qinv_s, double *s, int n) {
	double num, den;
	
	//num = exp(S[j + j * n] * (y_t)) * s[j];
	num = (S[j + j * n] * (y_t)) + log(s[j]);
	
	den = 0.0;
	Q += j;
	for(int i=0; i<n; i++) {
		den += *Q * exp(*(evals++) * (y_t)) * *(Qinv_s++);
		//den += exp(log(*Q) + (*(evals++) * (y_t)) + log(*(Qinv_s++))); // NO!!  Q might have negatives
		Q += n;
	}
	
	return(exp(num-log(den)));
	//return(num/den);
}


//// Sample a chain conditional on absorbing at an observation, returning only sufficient statistics of the chain
// y (input)
//     absorption time to be conditioned on
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
// Qinv_s (input)
//     precomputed n x 1 matrix Qinv %*% s
// n (input)
//     dimension
// PjQSjjLinv (input)
//     matrix where rows are P_{j\cdot} Q (\mathbf{I}(-S_{jj})+\Lambda)^{-1} for each j
// P (input)
//     n x n matrix of embedded transition probabilities, excluding absorbing moves
// res_z (output)
//     n dimensional array, filled with times spent in each state
// res_B (output)
//     integer, set to the state the sample chain started in
// res_N (output)
//     n x n matrix, filled with number of transitions between states (i != j).  The only non-zero diagonal was the exit-to-absorption state
// workD (output)
//     ?? element workspace
//
// return: sufficient statistics of the sample (z, B, N)
void LJMA_samplechain_Aslett2(double *y, double *pi, double *S, double *s, double *Q, double *evals, double *Qinv_s, int *n, double *PjQSjjLinv, double *evals_Sjj, double *P, double *res_z, int *res_B, int *res_N, double *workD) {
	LJMA_GetRNGstate();

	dotProd_pars pars;
	pars.PjQSjjLinv = PjQSjjLinv;
	pars.evals_Sjj = evals_Sjj;
	pars.diagConst = workD; workD += *n;
	pars.Ly_t = workD; workD += *n;
	pars.evals = evals;
	pars.Qinv_s = Qinv_s;
	pars.n = *n;
	
	// Initialise output vars
	*res_B = 0;
	int i;
	for(i=0; i<*n; i++) {
		res_z[i] = 0;
		for(int j=0; j<*n; j++) {
			res_N[i + j * *n] = 0;
		}
	}
	
	
	// Choose starting state ...
	double sofar = 0.0, target;
	target = runif(0.0, 1.0);
	*res_B = 0;
	while(sofar < target) {
		sofar += pi[(*res_B)++];
	}
	(*res_B)--;
	//Rprintf("%d - Start at: %d\n", *reverse, *res_B);
	
	
	// Run through the chain
	double t = 0.0, lastt = 0.0, d, y_t, Tol;
	int j = *res_B, lastj, Maxit;
	double *p;
	p = workD; workD += *n;
	while(TRUE) {
		y_t = *y-t;
		
		// Do we move from where we are, or is this the state from which we absorb?
		if(s[j] > 0.0) { // can't absorb from here if no exit allowed
			if(runif(0.0, 1.0) < LJMA_probAbsorb(y_t, j, S, Q, evals, Qinv_s, s, *n)) {
				break;
			}
		}
		
		lastt = t;
		lastj = j;
		
		/* D1 <- (1/(rep(-S[j,j],3)+L))*exp(L*(y-t))
		Dd <- (1/(rep(-S[j,j],3)+L))*exp(L*(y-t)-x*(L+rep(-S[j,j],3)))
		D2 <- (1/(rep(-S[j,j],3)+L))*exp(-rep((y-t)*(-S[j,j]),3))

		V <- Q %*% diag((D1*z-D2*z),3,3) %*% solve(Q)
		Ud <- Q %*% diag((D1-Dd),3,3) %*% solve(Q)

		P[j,]%*%Q%*%diag(D1*(z-1)-D2*z+Dd,3)%*%solve(Q)%*%s # It should be this
		 */		
		
		// Choose sojourn time
		pars.U = runif(0.0, 1.0);
		pars.j = j;
		pars.xi = -S[j + j * *n];
		for(i=0; i<*n; i++) {
			//pars.eLy_t[i] = exp(evals[i]*y_t); // Non-log scale
			pars.Ly_t[i] = evals[i]*y_t;
			//pars.diagConst[i] = (pars.U-1.0)*pars.eLy_t[i] - pars.U*exp(-y_t*pars.xi); // Non-log scale pre changing eLy_t
			//pars.diagConst[i] = (pars.U-1.0)*exp(pars.Ly_t[i]) - pars.U*exp(-y_t*pars.xi); // Non-log scale
			if( evals_Sjj[i + j * *n] != 0.0 ) {
				pars.diagConst[i] = -exp(log(1.0-pars.U)+pars.Ly_t[i]) - exp(log(pars.U)-y_t*pars.xi); // log scale
			} else {
				pars.diagConst[i] = -y_t*pars.U*exp(pars.Ly_t[i]);
			}
			//Rprintf("state %d, Ly_t = %e, diagConst = %e\n", i+1, pars.Ly_t[i], pars.diagConst[i]);
		}
		Tol = 0.0;
		Maxit = 10000;
		t += d = R_zeroin(0.0, y_t, &LJMA_dotProd, &pars, &Tol, &Maxit);
		//Rprintf("U=%lf, Tol=%e, Maxit=%d\n", pars.U, Tol, Maxit);
		
		// Make the state jump
		LJMA_moveMass(p, y_t-d, j, P, Q, evals, Qinv_s, *n, workD);

		target = runif(0.0, 1.0);
		sofar = 0.0;
		j = 0;
		while(sofar < target) {
			sofar += p[j++];
		}
		j--;
		
		//Rprintf("Sojourn: %lf, move %d -> %d\n", d, lastj, j);
		
		res_z[lastj] += d;
		res_N[lastj + j * *n]++;
	}
	
	//Rprintf("Exit from %d\n", j);
	
	res_N[j + j * *n]++;
	res_z[j] += *y-t;
	
	LJMA_GUI();
	LJMA_PutRNGstate();
}


//// (New Aslett METHOD) Sample a chain using directly conditional sampling, returning only sufficient statistics of the chain
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
// Qinv_s (input)
//     precomputed n x 1 matrix Qinv %*% s
// Qinv_1 (input)
//     precomputed n x 1 matrix Qinv %*% 1
// P (input)
//     n x n matrix of embedded transition probabilities, excluding absorbing moves
// Pfull (input)
//     n x (n+1) matrix of embedded chain transition probabilities, incluing option of absorption
// n (input)
//     dimension
// res_z (output)
//     n dimensional array, filled with sums of times spent in each state
// res_B (output)
//     n dimensional array, set to the number of times chains started in the given states
// res_N (output)
//     n x n matrix, filled with the sums of the number of transitions between states (i != j).  Diagonal entries are the exit-to-absorption counts
// workD (output)
//     ?? element workspace (SELF: )
// workI (output)
//     ?? element workspace (SELF: )
void LJMA_MHsample_Aslett2(double *y, int *censored, int *m, double *pi, double *S, double *s, double *Q, double *evals, double *Qinv_s, double *Qinv_1, double *P, double *Pfull, int *n, double *res_z, int *res_B, int *res_N, double *workD, int *workI) {
	// Initialise output vars
	int i,j,k,l;
	for(i=0; i<*n; i++) {
		res_B[i] = 0;
		res_z[i] = 0.0;
		for(j=0; j<*n; j++) {
			res_N[i + j * *n] = 0;
		}
	}
	
	// Precompute commonly used quantities
	char transN = 'N';
	double oneD = 1.0, zeroD = 0.0;
	
	double *evals_Sjj;
	evals_Sjj = workD; workD += *n * *n;
	for(j=0; j<*n; j++) {
		for(i=0; i<*n; i++) {
			//evals_Sjj[i + j * *n] = fabs((evals[i] - S[j + j * *n])/S[j + j * *n]) > DBL_EPSILON ? evals[i] - S[j + j * *n] : 0.0;
			evals_Sjj[i + j * *n] = fabs((evals[i] - S[j + j * *n])) > 3e-14 ? evals[i] - S[j + j * *n] : 0.0;
			//evals_Sjj[i + j * *n] = evals[i] - S[j + j * *n];
		}
	}
	
	double *PjQSjjLinv; // Used for jump time distribution
	PjQSjjLinv = workD; workD += *n * *n;
	F77_CALL(dgemm)(&transN, &transN, n, n, n, &oneD, P, n, Q, n, &zeroD, PjQSjjLinv, n); // First P_j Q ... but compute as (P_j Q)^T = Q^T P_j^T because want to store column wise for fast traversal
	for(j=0; j<*n; j++) { // ... then over -S_{jj} + \lambda_i
		for(i=0; i<*n; i++) {
			PjQSjjLinv[j + i * *n] = evals_Sjj[i + j * *n] != 0.0 ? PjQSjjLinv[j + i * *n]/evals_Sjj[i + j * *n] : PjQSjjLinv[j + i * *n];
			//Rprintf("evals_Sjj[%d,%d] = %e\nPjQSjjLinv[%d,%d] = %e\n",i,j,evals_Sjj[i + j * *n],j,i,PjQSjjLinv[j + i * *n]);//, -S[j + j * *n], (-S[j + j * *n] + evals[i]), fabs((-S[j + j * *n] + evals[i])/S[j + j * *n]) > 2*DBL_EPSILON);
		}
	}
	
	
	// Create memory for generated sufficient stats
	double *z;
	int B, *N;
	z = workD; workD += *n;
	N = workI; workI += *n * *n;
	
	
	// Off we go
	double *y_p;
	int *censored_p, reverse=0, res_pre;
	y_p = y;
	censored_p = censored;
	for(i=0; i<*m; i++) {
		if(*censored_p) {
			LJMA_samplechain(y_p, censored_p, pi, S, Q, evals, Qinv_1, P, Pfull, n, &reverse, NULL, z, &B, N, &res_pre, workD);
		} else {
			LJMA_samplechain_Aslett2(y_p, pi, S, s, Q, evals, Qinv_s, n, PjQSjjLinv, evals_Sjj, P, z, &B, N, workD);
		}
		
		// Now add this step to the running total
		res_B[B]++;
		for(k=0; k<*n; k++) {
			res_z[k] += z[k];
			for(l=0; l<*n; l++) {
				res_N[k + l * *n] += N[k + l * *n];
			}
		}
		
		y_p++;
		censored_p++;
	}
}


