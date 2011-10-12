/*
 *  Simulate_AbsCTMC_eq_Bladt_MHRS.h
 *  
 *
 *  Created by Louis Aslett (www.louisaslett.com)
 *  
 *  Functions to simulate from absorbing Continuous-time Markov Chain conditional on a given absorbtion time.
 *
 */

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
void LJMA_MHsample_Bladt(double *y, int *censored, int *m, double *pi, double *S, double *s, double *Pfull, int *n, int *iter, double *res_z, int *res_B, int *res_N, double *workD, int *workI);
