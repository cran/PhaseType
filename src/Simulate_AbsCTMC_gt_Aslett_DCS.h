/*
 *  Simulate_AbsCTMC_gt_Aslett_DCS.h
 *  
 *
 *  Created by Louis Aslett (www.louisaslett.com)
 *  
 *  Functions to simulate from absorbing Continuous-time Markov Chain conditional on exceeding a given absorbtion time.
 *
 */

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
//     5n element workspace
//
// return: sufficient statistics of the sample (z, B, N)
void LJMA_samplechain(double *y, int *censored, double *pi, double *S, double *Q, double *evals, double *Qinv_1, double *P, double *Pfull, int *n, int *reverse, double *piR, double *res_z, int *res_B, int *res_N, int *res_pre, double *workD);
