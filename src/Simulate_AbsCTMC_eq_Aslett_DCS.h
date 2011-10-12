/*
 *  Simulate_AbsCTMC_eq_Aslett_DCS.h
 *  
 *
 *  Created by Louis Aslett (www.louisaslett.com)
 *  
 *  Functions to simulate from absorbing Continuous-time Markov Chain conditional on a given absorbtion time.
 *
 */

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
void LJMA_MHsample_Aslett(double *y, int *censored, int *m, double *pi, double *S, double *s, double *P, double *Pfull, double *Q, double *evals, double *Qinv_1, double *Qinv, int *n, int *iter, int *reverse, double *res_z, int *res_B, int *res_N, double *workD, int *workI);

