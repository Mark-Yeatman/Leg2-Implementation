/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: SLIP_KPBC_v4.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 26-Sep-2019 12:16:46
 */

#ifndef SLIP_KPBC_V4_H
#define SLIP_KPBC_V4_H

/* Include Files */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "SLIP_KPBC_v4_types.h"

/* Function Declarations */
extern void SLIP_KPBC_v4(float IMU_pitch, float Knee_motor_position, float
  Knee_joint_position, float Ankle_motor_position, float Ankle_joint_position,
  float Iteration, float Iteration_time, float Knee_torque_sensor, float
  Ankle_torque_sensor, float Load_cell_x_force, float Load_cell_y_force, float
  Load_cell_z_force, float Load_cell_x_moment, float Load_cell_y_moment, float
  Load_cell_z_moment, float kp_knee, float kd_knee, float kp_ankle, float
  kd_ankle, float ankle_ind_of_hip, float knee_ind_of_hip, float ankle_des_in,
  float knee_des_in, float time_in, float filter_coeff, float IMU_filter_coeff,
  float q_h_0, float q_h_min, float c, float s_po, float FC, float lf, float la,
  float ls, float lt, float k, float L0, float *Knee_torque_command, float
  *Ankle_torque_command, float *deltaL, float *hip_pos, float *t_out, float *dt);
extern void knee_pos_prev_not_empty_init(void);
extern void t_prev_not_empty_init(void);

#endif

/*
 * File trailer for SLIP_KPBC_v4.h
 *
 * [EOF]
 */
