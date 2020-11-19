/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_SLIP_KPBC_api.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 31-Aug-2020 12:07:22
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_SLIP_KPBC_api.h"
#include "_coder_SLIP_KPBC_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131435U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "SLIP_KPBC",                         /* fFunctionName */
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
 * Arguments    : const mxArray * const prhs[40]
 *                const mxArray *plhs[27]
 * Return Type  : void
 */
void SLIP_KPBC_api(const mxArray * const prhs[40], const mxArray *plhs[27])
{
  real32_T IMU_pitch;
  real32_T Knee_joint_position;
  real32_T Ankle_joint_position;
  real32_T time_in;
  real32_T dt;
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
  real32_T ankle_des_in;
  real32_T knee_des_in;
  real32_T vel_filter_coeff;
  real32_T KPBC_filter_coeff;
  real32_T SLIP_ON;
  real32_T lt;
  real32_T k;
  real32_T d;
  real32_T L0;
  real32_T KPBC_ON;
  real32_T KPBC_max_torque_rate;
  real32_T pbc_gain_knee;
  real32_T md;
  real32_T Eref;
  real32_T knee_stop_low;
  real32_T knee_stop_high;
  real32_T ankle_stop_low;
  real32_T ankle_stop_high;
  real32_T max_torque;
  real32_T v0;
  real32_T Fric_Comp;
  real32_T F_thresh;
  real32_T Command_State;
  real32_T KPBC_max_torque;
  real32_T Joint_Bio_Sat;
  real32_T Knee_torque_command;
  real32_T Ankle_torque_command;
  real32_T deltaL;
  real32_T hip_pos;
  real32_T Esys;
  real32_T Esys_integrate_out;
  real32_T U_S_KNEE;
  real32_T U_S_ANKLE;
  real32_T COPFX;
  real32_T U_LIN_DAMP_A;
  real32_T U_STOP_K;
  real32_T U_STOP_A;
  real32_T U_PBC_K;
  real32_T U_PBC_A;
  real32_T knee_des_out;
  real32_T ankle_des_out;
  real32_T foot_contact;
  real32_T stance;
  real32_T swing;
  real32_T phase_var_out;
  real32_T IMU_LIVE_OUT;
  real32_T StanceGain;
  real32_T SwingGain;
  real32_T knee_joint_vel;
  real32_T ankle_joint_vel;
  real32_T hip_vel;
  real32_T PushOffOut;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  IMU_pitch = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "IMU_pitch");
  Knee_joint_position = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]),
    "Knee_joint_position");
  Ankle_joint_position = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]),
    "Ankle_joint_position");
  time_in = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "time_in");
  dt = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "dt");
  Load_cell_x_force = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]),
    "Load_cell_x_force");
  Load_cell_y_force = emlrt_marshallIn(&st, emlrtAliasP(prhs[6]),
    "Load_cell_y_force");
  Load_cell_z_force = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]),
    "Load_cell_z_force");
  Load_cell_x_moment = emlrt_marshallIn(&st, emlrtAliasP(prhs[8]),
    "Load_cell_x_moment");
  Load_cell_y_moment = emlrt_marshallIn(&st, emlrtAliasP(prhs[9]),
    "Load_cell_y_moment");
  Load_cell_z_moment = emlrt_marshallIn(&st, emlrtAliasP(prhs[10]),
    "Load_cell_z_moment");
  kp_knee = emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "kp_knee");
  kd_knee = emlrt_marshallIn(&st, emlrtAliasP(prhs[12]), "kd_knee");
  kp_ankle = emlrt_marshallIn(&st, emlrtAliasP(prhs[13]), "kp_ankle");
  kd_ankle = emlrt_marshallIn(&st, emlrtAliasP(prhs[14]), "kd_ankle");
  ankle_des_in = emlrt_marshallIn(&st, emlrtAliasP(prhs[15]), "ankle_des_in");
  knee_des_in = emlrt_marshallIn(&st, emlrtAliasP(prhs[16]), "knee_des_in");
  vel_filter_coeff = emlrt_marshallIn(&st, emlrtAliasP(prhs[17]),
    "vel_filter_coeff");
  KPBC_filter_coeff = emlrt_marshallIn(&st, emlrtAliasP(prhs[18]),
    "KPBC_filter_coeff");
  SLIP_ON = emlrt_marshallIn(&st, emlrtAliasP(prhs[19]), "SLIP_ON");
  lt = emlrt_marshallIn(&st, emlrtAliasP(prhs[20]), "lt");
  k = emlrt_marshallIn(&st, emlrtAliasP(prhs[21]), "k");
  d = emlrt_marshallIn(&st, emlrtAliasP(prhs[22]), "d");
  L0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[23]), "L0");
  KPBC_ON = emlrt_marshallIn(&st, emlrtAliasP(prhs[24]), "KPBC_ON");
  KPBC_max_torque_rate = emlrt_marshallIn(&st, emlrtAliasP(prhs[25]),
    "KPBC_max_torque_rate");
  pbc_gain_knee = emlrt_marshallIn(&st, emlrtAliasP(prhs[26]), "pbc_gain_knee");
  md = emlrt_marshallIn(&st, emlrtAliasP(prhs[27]), "md");
  Eref = emlrt_marshallIn(&st, emlrtAliasP(prhs[28]), "Eref");
  knee_stop_low = emlrt_marshallIn(&st, emlrtAliasP(prhs[29]), "knee_stop_low");
  knee_stop_high = emlrt_marshallIn(&st, emlrtAliasP(prhs[30]), "knee_stop_high");
  ankle_stop_low = emlrt_marshallIn(&st, emlrtAliasP(prhs[31]), "ankle_stop_low");
  ankle_stop_high = emlrt_marshallIn(&st, emlrtAliasP(prhs[32]),
    "ankle_stop_high");
  max_torque = emlrt_marshallIn(&st, emlrtAliasP(prhs[33]), "max_torque");
  v0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[34]), "v0");
  Fric_Comp = emlrt_marshallIn(&st, emlrtAliasP(prhs[35]), "Fric_Comp");
  F_thresh = emlrt_marshallIn(&st, emlrtAliasP(prhs[36]), "F_thresh");
  Command_State = emlrt_marshallIn(&st, emlrtAliasP(prhs[37]), "Command_State");
  KPBC_max_torque = emlrt_marshallIn(&st, emlrtAliasP(prhs[38]),
    "KPBC_max_torque");
  Joint_Bio_Sat = emlrt_marshallIn(&st, emlrtAliasP(prhs[39]), "Joint_Bio_Sat");

  /* Invoke the target function */
  SLIP_KPBC(IMU_pitch, Knee_joint_position, Ankle_joint_position, time_in, dt,
            Load_cell_x_force, Load_cell_y_force, Load_cell_z_force,
            Load_cell_x_moment, Load_cell_y_moment, Load_cell_z_moment, kp_knee,
            kd_knee, kp_ankle, kd_ankle, ankle_des_in, knee_des_in,
            vel_filter_coeff, KPBC_filter_coeff, SLIP_ON, lt, k, d, L0, KPBC_ON,
            KPBC_max_torque_rate, pbc_gain_knee, md, Eref, knee_stop_low,
            knee_stop_high, ankle_stop_low, ankle_stop_high, max_torque, v0,
            Fric_Comp, F_thresh, Command_State, KPBC_max_torque, Joint_Bio_Sat,
            &Knee_torque_command, &Ankle_torque_command, &deltaL, &hip_pos,
            &Esys, &Esys_integrate_out, &U_S_KNEE, &U_S_ANKLE, &COPFX,
            &U_LIN_DAMP_A, &U_STOP_K, &U_STOP_A, &U_PBC_K, &U_PBC_A,
            &knee_des_out, &ankle_des_out, &foot_contact, &stance, &swing,
            &phase_var_out, &IMU_LIVE_OUT, &StanceGain, &SwingGain,
            &knee_joint_vel, &ankle_joint_vel, &hip_vel, &PushOffOut);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(Knee_torque_command);
  plhs[1] = emlrt_marshallOut(Ankle_torque_command);
  plhs[2] = emlrt_marshallOut(deltaL);
  plhs[3] = emlrt_marshallOut(hip_pos);
  plhs[4] = emlrt_marshallOut(Esys);
  plhs[5] = emlrt_marshallOut(Esys_integrate_out);
  plhs[6] = emlrt_marshallOut(U_S_KNEE);
  plhs[7] = emlrt_marshallOut(U_S_ANKLE);
  plhs[8] = emlrt_marshallOut(COPFX);
  plhs[9] = emlrt_marshallOut(U_LIN_DAMP_A);
  plhs[10] = emlrt_marshallOut(U_STOP_K);
  plhs[11] = emlrt_marshallOut(U_STOP_A);
  plhs[12] = emlrt_marshallOut(U_PBC_K);
  plhs[13] = emlrt_marshallOut(U_PBC_A);
  plhs[14] = emlrt_marshallOut(knee_des_out);
  plhs[15] = emlrt_marshallOut(ankle_des_out);
  plhs[16] = emlrt_marshallOut(foot_contact);
  plhs[17] = emlrt_marshallOut(stance);
  plhs[18] = emlrt_marshallOut(swing);
  plhs[19] = emlrt_marshallOut(phase_var_out);
  plhs[20] = emlrt_marshallOut(IMU_LIVE_OUT);
  plhs[21] = emlrt_marshallOut(StanceGain);
  plhs[22] = emlrt_marshallOut(SwingGain);
  plhs[23] = emlrt_marshallOut(knee_joint_vel);
  plhs[24] = emlrt_marshallOut(ankle_joint_vel);
  plhs[25] = emlrt_marshallOut(hip_vel);
  plhs[26] = emlrt_marshallOut(PushOffOut);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void SLIP_KPBC_atexit(void)
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
  SLIP_KPBC_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void SLIP_KPBC_initialize(void)
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
void SLIP_KPBC_terminate(void)
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
 * File trailer for _coder_SLIP_KPBC_api.c
 *
 * [EOF]
 */
