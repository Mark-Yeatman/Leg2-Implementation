/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_SLIP_api.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 15-Oct-2019 15:51:04
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_SLIP_api.h"
#include "_coder_SLIP_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131435U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "SLIP",                              /* fFunctionName */
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
static real32_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *IMU_pitch,
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
 *                const mxArray *IMU_pitch
 *                const char_T *identifier
 * Return Type  : real32_T
 */
static real32_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *IMU_pitch,
  const char_T *identifier)
{
  real32_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(IMU_pitch), &thisId);
  emlrtDestroyArray(&IMU_pitch);
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
 * Arguments    : const mxArray * const prhs[39]
 *                const mxArray *plhs[6]
 * Return Type  : void
 */
void SLIP_api(const mxArray * const prhs[39], const mxArray *plhs[6])
{
  real32_T IMU_pitch;
  real32_T Knee_motor_position;
  real32_T Knee_joint_position;
  real32_T Ankle_motor_position;
  real32_T Ankle_joint_position;
  real32_T Iteration;
  real32_T Iteration_time;
  real32_T Knee_torque_sensor;
  real32_T Ankle_torque_sensor;
  real32_T Load_cell_x_force;
  real32_T Load_cell_y_force;
  real32_T Load_cell_z_force;
  real32_T Load_cell_x_moment;
  real32_T Load_cell_y_moment;
  real32_T Load_cell_z_moment;
  real32_T kp_knee;
  real32_T kd_knee;
  real32_T kp_ankle;
  real32_T kd_ankle;
  real32_T ankle_ind_of_hip;
  real32_T knee_ind_of_hip;
  real32_T ankle_des_in;
  real32_T knee_des_in;
  real32_T time_in;
  real32_T filter_coeff;
  real32_T IMU_filter_coeff;
  real32_T q_h_0;
  real32_T q_h_min;
  real32_T c;
  real32_T s_po;
  real32_T FC;
  real32_T lf;
  real32_T la;
  real32_T ls;
  real32_T lt;
  real32_T k;
  real32_T d;
  real32_T L0;
  real32_T SLIP_ON;
  real32_T Knee_torque_command;
  real32_T Ankle_torque_command;
  real32_T deltaL;
  real32_T hip_pos;
  real32_T t_out;
  real32_T dt;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  IMU_pitch = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "IMU_pitch");
  Knee_motor_position = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]),
    "Knee_motor_position");
  Knee_joint_position = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]),
    "Knee_joint_position");
  Ankle_motor_position = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]),
    "Ankle_motor_position");
  Ankle_joint_position = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]),
    "Ankle_joint_position");
  Iteration = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "Iteration");
  Iteration_time = emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "Iteration_time");
  Knee_torque_sensor = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]),
    "Knee_torque_sensor");
  Ankle_torque_sensor = emlrt_marshallIn(&st, emlrtAliasP(prhs[8]),
    "Ankle_torque_sensor");
  Load_cell_x_force = emlrt_marshallIn(&st, emlrtAliasP(prhs[9]),
    "Load_cell_x_force");
  Load_cell_y_force = emlrt_marshallIn(&st, emlrtAliasP(prhs[10]),
    "Load_cell_y_force");
  Load_cell_z_force = emlrt_marshallIn(&st, emlrtAliasP(prhs[11]),
    "Load_cell_z_force");
  Load_cell_x_moment = emlrt_marshallIn(&st, emlrtAliasP(prhs[12]),
    "Load_cell_x_moment");
  Load_cell_y_moment = emlrt_marshallIn(&st, emlrtAliasP(prhs[13]),
    "Load_cell_y_moment");
  Load_cell_z_moment = emlrt_marshallIn(&st, emlrtAliasP(prhs[14]),
    "Load_cell_z_moment");
  kp_knee = emlrt_marshallIn(&st, emlrtAliasP(prhs[15]), "kp_knee");
  kd_knee = emlrt_marshallIn(&st, emlrtAliasP(prhs[16]), "kd_knee");
  kp_ankle = emlrt_marshallIn(&st, emlrtAliasP(prhs[17]), "kp_ankle");
  kd_ankle = emlrt_marshallIn(&st, emlrtAliasP(prhs[18]), "kd_ankle");
  ankle_ind_of_hip = emlrt_marshallIn(&st, emlrtAliasP(prhs[19]),
    "ankle_ind_of_hip");
  knee_ind_of_hip = emlrt_marshallIn(&st, emlrtAliasP(prhs[20]),
    "knee_ind_of_hip");
  ankle_des_in = emlrt_marshallIn(&st, emlrtAliasP(prhs[21]), "ankle_des_in");
  knee_des_in = emlrt_marshallIn(&st, emlrtAliasP(prhs[22]), "knee_des_in");
  time_in = emlrt_marshallIn(&st, emlrtAliasP(prhs[23]), "time_in");
  filter_coeff = emlrt_marshallIn(&st, emlrtAliasP(prhs[24]), "filter_coeff");
  IMU_filter_coeff = emlrt_marshallIn(&st, emlrtAliasP(prhs[25]),
    "IMU_filter_coeff");
  q_h_0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[26]), "q_h_0");
  q_h_min = emlrt_marshallIn(&st, emlrtAliasP(prhs[27]), "q_h_min");
  c = emlrt_marshallIn(&st, emlrtAliasP(prhs[28]), "c");
  s_po = emlrt_marshallIn(&st, emlrtAliasP(prhs[29]), "s_po");
  FC = emlrt_marshallIn(&st, emlrtAliasP(prhs[30]), "FC");
  lf = emlrt_marshallIn(&st, emlrtAliasP(prhs[31]), "lf");
  la = emlrt_marshallIn(&st, emlrtAliasP(prhs[32]), "la");
  ls = emlrt_marshallIn(&st, emlrtAliasP(prhs[33]), "ls");
  lt = emlrt_marshallIn(&st, emlrtAliasP(prhs[34]), "lt");
  k = emlrt_marshallIn(&st, emlrtAliasP(prhs[35]), "k");
  d = emlrt_marshallIn(&st, emlrtAliasP(prhs[36]), "d");
  L0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[37]), "L0");
  SLIP_ON = emlrt_marshallIn(&st, emlrtAliasP(prhs[38]), "SLIP_ON");

  /* Invoke the target function */
  SLIP(IMU_pitch, Knee_motor_position, Knee_joint_position, Ankle_motor_position,
       Ankle_joint_position, Iteration, Iteration_time, Knee_torque_sensor,
       Ankle_torque_sensor, Load_cell_x_force, Load_cell_y_force,
       Load_cell_z_force, Load_cell_x_moment, Load_cell_y_moment,
       Load_cell_z_moment, kp_knee, kd_knee, kp_ankle, kd_ankle,
       ankle_ind_of_hip, knee_ind_of_hip, ankle_des_in, knee_des_in, time_in,
       filter_coeff, IMU_filter_coeff, q_h_0, q_h_min, c, s_po, FC, lf, la, ls,
       lt, k, d, L0, SLIP_ON, &Knee_torque_command, &Ankle_torque_command,
       &deltaL, &hip_pos, &t_out, &dt);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(Knee_torque_command);
  plhs[1] = emlrt_marshallOut(Ankle_torque_command);
  plhs[2] = emlrt_marshallOut(deltaL);
  plhs[3] = emlrt_marshallOut(hip_pos);
  plhs[4] = emlrt_marshallOut(t_out);
  plhs[5] = emlrt_marshallOut(dt);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void SLIP_atexit(void)
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
  SLIP_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void SLIP_initialize(void)
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
void SLIP_terminate(void)
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
 * File trailer for _coder_SLIP_api.c
 *
 * [EOF]
 */
