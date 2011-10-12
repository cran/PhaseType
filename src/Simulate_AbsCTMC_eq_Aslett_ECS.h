/*
 *  Simulate_AbsCTMC_eq_Aslett_ECS.h
 *  
 *
 *  Created by Louis Aslett (www.louisaslett.com)
 *  
 *  Functions to simulate from absorbing Continuous-time Markov Chain conditional on a given absorbtion time.
 *
 */

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
void LJMA_MHsample_Aslett2(double *y, int *censored, int *m, double *pi, double *S, double *s, double *Q, double *evals, double *Qinv_s, double *Qinv_1, double *P, double *Pfull, int *n, double *res_z, int *res_B, int *res_N, double *workD, int *workI);
	
