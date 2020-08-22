/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_SLIP_api.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 15-Oct-2019 15:51:04
 */

#ifndef _CODER_SLIP_API_H
#define _CODER_SLIP_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_SLIP_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void SLIP(real32_T IMU_pitch, real32_T Knee_motor_position, real32_T
                 Knee_joint_position, real32_T Ankle_motor_position, real32_T
                 Ankle_joint_position, real32_T Iteration, real32_T
                 Iteration_time, real32_T Knee_torque_sensor, real32_T
                 Ankle_torque_sensor, real32_T Load_cell_x_force, real32_T
                 Load_cell_y_force, real32_T Load_cell_z_force, real32_T
                 Load_cell_x_moment, real32_T Load_cell_y_moment, real32_T
                 Load_cell_z_moment, real32_T kp_knee, real32_T kd_knee,
                 real32_T kp_ankle, real32_T kd_ankle, real32_T ankle_ind_of_hip,
                 real32_T knee_ind_of_hip, real32_T ankle_des_in, real32_T
                 knee_des_in, real32_T time_in, real32_T filter_coeff, real32_T
                 IMU_filter_coeff, real32_T q_h_0, real32_T q_h_min, real32_T c,
                 real32_T s_po, real32_T FC, real32_T lf, real32_T la, real32_T
                 ls, real32_T lt, real32_T k, real32_T d, real32_T L0, real32_T
                 SLIP_ON, real32_T *Knee_torque_command, real32_T
                 *Ankle_torque_command, real32_T *deltaL, real32_T *hip_pos,
                 real32_T *t_out, real32_T *dt);
extern void SLIP_api(const mxArray * const prhs[39], const mxArray *plhs[6]);
extern void SLIP_atexit(void);
extern void SLIP_initialize(void);
extern void SLIP_terminate(void);
extern void SLIP_xil_terminate(void);

#endif

/*
 * File trailer for _coder_SLIP_api.h
 *
 * [EOF]
 */
