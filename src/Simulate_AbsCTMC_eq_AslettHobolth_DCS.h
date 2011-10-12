/*
 *  Simulate_AbsCTMC_eq_AslettHobolth_DCS.h
 *  
 *
 *  Created by Louis Aslett (www.louisaslett.com)
 *  
 *  Functions to simulate from absorbing Continuous-time Markov Chain conditional on a given absorbtion time.
 *
 */

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
// bvec (input)
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
void LJMA_MHsample_Hobolth2(double *y, int *censored, int *m, double *pi, double *S, double *s, double *Q, double *evals, double *Qinv_b, double *bvec, double *Qinv, int *n, int *iter, double *res_z, int *res_B, int *res_N, double *workD, int *workI);
	
