#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "PHT_MCMC_Aslett.h"

static R_NativePrimitiveArgType C_type1[] = {
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP
};

static const R_CMethodDef cMethods[] = {
  {"LJMA_Gibbs", (DL_FUNC) &LJMA_Gibbs, 15, C_type1},
  {NULL, NULL, 0, NULL}
};

void R_init_PhaseType(DllInfo *info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
