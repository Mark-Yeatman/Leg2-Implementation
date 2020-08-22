/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_SLIP_KPBC_v4_mex.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 26-Sep-2019 12:16:46
 */

/* Include Files */
#include "_coder_SLIP_KPBC_v4_api.h"
#include "_coder_SLIP_KPBC_v4_mex.h"

/* Function Declarations */
static void SLIP_KPBC_v4_mexFunction(int32_T nlhs, mxArray *plhs[6], int32_T
  nrhs, const mxArray *prhs[37]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                const mxArray *plhs[6]
 *                int32_T nrhs
 *                const mxArray *prhs[37]
 * Return Type  : void
 */
static void SLIP_KPBC_v4_mexFunction(int32_T nlhs, mxArray *plhs[6], int32_T
  nrhs, const mxArray *prhs[37])
{
  int32_T n;
  const mxArray *inputs[37];
  const mxArray *outputs[6];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 37) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 37, 4,
                        12, "SLIP_KPBC_v4");
  }

  if (nlhs > 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 12,
                        "SLIP_KPBC_v4");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
  }

  /* Call the function. */
  SLIP_KPBC_v4_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  SLIP_KPBC_v4_terminate();
}

/*
 * Arguments    : int32_T nlhs
 *                const mxArray * const plhs[]
 *                int32_T nrhs
 *                const mxArray * const prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(SLIP_KPBC_v4_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  SLIP_KPBC_v4_initialize();

  /* Dispatch the entry-point. */
  SLIP_KPBC_v4_mexFunction(nlhs, plhs, nrhs, prhs);
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_SLIP_KPBC_v4_mex.c
 *
 * [EOF]
 */
