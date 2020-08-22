/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_PDControlTest_api.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 03-Oct-2019 15:42:50
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_PDControlTest_api.h"
#include "_coder_PDControlTest_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131435U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "PDControlTest",                     /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real32_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real32_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);
static real32_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *knee_pos,
  const char_T *identifier);
static const mxArray *emlrt_marshallOut(const real32_T u);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real32_T
 */
static real32_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real32_T y;
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real32_T
 */
static real32_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  real32_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "single|double", false, 0U, &dims);
  emlrtImportArrayR2015b(sp, src, &ret, 4, false);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *knee_pos
 *                const char_T *identifier
 * Return Type  : real32_T
 */
static real32_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *knee_pos,
  const char_T *identifier)
{
  real32_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(knee_pos), &thisId);
  emlrtDestroyArray(&knee_pos);
  return y;
}

/*
 * Arguments    : const real32_T u
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real32_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
  *(real32_T *)mxGetData(m0) = u;
  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const mxArray * const prhs[11]
 *                const mxArray *plhs[4]
 * Return Type  : void
 */
void PDControlTest_api(const mxArray * const prhs[11], const mxArray *plhs[4])
{
  real32_T knee_pos;
  real32_T ankle_pos;
  real32_T knee_des;
  real32_T ankle_des;
  real32_T dt;
  real32_T kp_k;
  real32_T kd_k;
  real32_T kp_a;
  real32_T kd_a;
  real32_T time_in;
  real32_T filter_coeff;
  real32_T u_k;
  real32_T u_a;
  real32_T t_out;
  real32_T dt_out;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  knee_pos = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "knee_pos");
  ankle_pos = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "ankle_pos");
  knee_des = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "knee_des");
  ankle_des = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "ankle_des");
  dt = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "dt");
  kp_k = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "kp_k");
  kd_k = emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "kd_k");
  kp_a = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "kp_a");
  kd_a = emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "kd_a");
  time_in = emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "time_in");
  filter_coeff = emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "filter_coeff");

  /* Invoke the target function */
  PDControlTest(knee_pos, ankle_pos, knee_des, ankle_des, dt, kp_k, kd_k, kp_a,
                kd_a, time_in, filter_coeff, &u_k, &u_a, &t_out, &dt_out);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(u_k);
  plhs[1] = emlrt_marshallOut(u_a);
  plhs[2] = emlrt_marshallOut(t_out);
  plhs[3] = emlrt_marshallOut(dt_out);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void PDControlTest_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  PDControlTest_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void PDControlTest_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void PDControlTest_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_PDControlTest_api.c
 *
 * [EOF]
 */
