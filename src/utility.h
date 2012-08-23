#include <R.h>
#include <R_ext/Utils.h>

extern double *LJMA_LAPACK_work;
extern int LJMA_LAPACK_lwork;

extern int LJMA_counter;

// The next two functions are for allocating the LAPACK work space when calling certain functions directly from R instead of via the outer Gibbs call
void LJMA_LAPACKspace(int *n);
void LJMA_LAPACKspaceFree();

// Then we can get at the performance counter here
void LJMA_setCounter(int *counter);
void LJMA_getCounter(int *counter);

// Functions to track random number generator state
extern int LJMA_RNG;
/*R_INLINE void LJMA_GetRNGstate() {
	if(LJMA_RNG++ == 0) GetRNGstate();
}
R_INLINE void LJMA_PutRNGstate() {
	if(--LJMA_RNG == 0) PutRNGstate();
}
R_INLINE void LJMA_GUI() {
	R_FlushConsole();
	#if ( defined(HAVE_AQUA) || defined(Win32) )
	R_CheckUserInterrupt();
	#endif
}*/
#define LJMA_GetRNGstate() { if(LJMA_RNG++ == 0) GetRNGstate(); }
#define LJMA_PutRNGstate() { if(--LJMA_RNG == 0) PutRNGstate(); }
#if ( defined(HAVE_AQUA) || defined(Win32) )
	#define LJMA_GUI() { \
		R_FlushConsole(); \
		R_CheckUserInterrupt(); \
	}
#else
	#define LJMA_GUI() { \
		R_FlushConsole(); \
	}
#endif

int LJMA_inverse(double *A, int *n, int *workI);

int LJMA_eigen(int *n, double *S, double *evals, double *Q, double *Qinv, double *workD, int *workI);

double Find0(double ax, double bx, double (*f)(double, void *), void *info, double *Tol, int *Maxit);
double Find02(double ax, double bx, double fa, double fb, double (*f)(double, void *), void *info, double *Tol, int *Maxit);
