/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_PDControlTest_api.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 03-Oct-2019 15:42:50
 */

#ifndef _CODER_PDCONTROLTEST_API_H
#define _CODER_PDCONTROLTEST_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_PDControlTest_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void PDControlTest(real32_T knee_pos, real32_T ankle_pos, real32_T
  knee_des, real32_T ankle_des, real32_T dt, real32_T kp_k, real32_T kd_k,
  real32_T kp_a, real32_T kd_a, real32_T time_in, real32_T filter_coeff,
  real32_T *u_k, real32_T *u_a, real32_T *t_out, real32_T *dt_out);
extern void PDControlTest_api(const mxArray * const prhs[11], const mxArray
  *plhs[4]);
extern void PDControlTest_atexit(void);
extern void PDControlTest_initialize(void);
extern void PDControlTest_terminate(void);
extern void PDControlTest_xil_terminate(void);

#endif

/*
 * File trailer for _coder_PDControlTest_api.h
 *
 * [EOF]
 */
