/*
 *  Simulate_AbsCTMC_gt_Bladt_MHRS.h
 *  
 *
 *  Created by Louis Aslett.
 *  
 *  Functions to simulate from absorbing Continuous-time Markov Chain conditional on exceeding a given absorbtion time.
 *
 */

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
void LJMA_samplechain_Bladt(double *y, int *censored, double *pi, double *S, double *Pfull, int *n, double *absorbed, double *res_z, int *res_B, int *res_N, int *res_pre, double *workD, int *workI);
