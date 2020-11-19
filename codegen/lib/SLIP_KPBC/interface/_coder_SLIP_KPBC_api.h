/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_SLIP_KPBC_api.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 31-Aug-2020 12:07:22
 */

#ifndef _CODER_SLIP_KPBC_API_H
#define _CODER_SLIP_KPBC_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_SLIP_KPBC_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void SLIP_KPBC(real32_T IMU_pitch, real32_T Knee_joint_position, real32_T
                      Ankle_joint_position, real32_T time_in, real32_T dt,
                      real32_T Load_cell_x_force, real32_T Load_cell_y_force,
                      real32_T Load_cell_z_force, real32_T Load_cell_x_moment,
                      real32_T Load_cell_y_moment, real32_T Load_cell_z_moment,
                      real32_T kp_knee, real32_T kd_knee, real32_T kp_ankle,
                      real32_T kd_ankle, real32_T ankle_des_in, real32_T
                      knee_des_in, real32_T vel_filter_coeff, real32_T
                      KPBC_filter_coeff, real32_T SLIP_ON, real32_T lt, real32_T
                      k, real32_T d, real32_T L0, real32_T KPBC_ON, real32_T
                      KPBC_max_torque_rate, real32_T pbc_gain_knee, real32_T md,
                      real32_T Eref, real32_T knee_stop_low, real32_T
                      knee_stop_high, real32_T ankle_stop_low, real32_T
                      ankle_stop_high, real32_T max_torque, real32_T v0,
                      real32_T Fric_Comp, real32_T F_thresh, real32_T
                      Command_State, real32_T KPBC_max_torque, real32_T
                      Joint_Bio_Sat, real32_T *Knee_torque_command, real32_T
                      *Ankle_torque_command, real32_T *deltaL, real32_T *hip_pos,
                      real32_T *Esys, real32_T *Esys_integrate_out, real32_T
                      *U_S_KNEE, real32_T *U_S_ANKLE, real32_T *COPFX, real32_T *
                      U_LIN_DAMP_A, real32_T *U_STOP_K, real32_T *U_STOP_A,
                      real32_T *U_PBC_K, real32_T *U_PBC_A, real32_T
                      *knee_des_out, real32_T *ankle_des_out, real32_T
                      *foot_contact, real32_T *stance, real32_T *swing, real32_T
                      *phase_var_out, real32_T *IMU_LIVE_OUT, real32_T
                      *StanceGain, real32_T *SwingGain, real32_T *knee_joint_vel,
                      real32_T *ankle_joint_vel, real32_T *hip_vel, real32_T
                      *PushOffOut);
extern void SLIP_KPBC_api(const mxArray * const prhs[40], const mxArray *plhs[27]);
extern void SLIP_KPBC_atexit(void);
extern void SLIP_KPBC_initialize(void);
extern void SLIP_KPBC_terminate(void);
extern void SLIP_KPBC_xil_terminate(void);

#endif

/*
 * File trailer for _coder_SLIP_KPBC_api.h
 *
 * [EOF]
 */
