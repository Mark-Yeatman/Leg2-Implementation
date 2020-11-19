/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: SLIP_KPBC.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 31-Aug-2020 12:07:22
 */

#ifndef SLIP_KPBC_H
#define SLIP_KPBC_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "SLIP_KPBC_types.h"

/* Function Declarations */
extern void SLIP_KPBC(float IMU_pitch, float Knee_joint_position, float
                      Ankle_joint_position, float time_in, float dt, float
                      Load_cell_x_force, float Load_cell_y_force, float
                      Load_cell_z_force, float Load_cell_x_moment, float
                      Load_cell_y_moment, float Load_cell_z_moment, float
                      kp_knee, float kd_knee, float kp_ankle, float kd_ankle,
                      float ankle_des_in, float knee_des_in, float
                      vel_filter_coeff, float KPBC_filter_coeff, float SLIP_ON,
                      float lt, float k, float d, float L0, float KPBC_ON, float
                      KPBC_max_torque_rate, float pbc_gain_knee, float md, float
                      Eref, float knee_stop_low, float knee_stop_high, float
                      ankle_stop_low, float ankle_stop_high, float max_torque,
                      float v0, float Fric_Comp, float F_thresh, float
                      Command_State, float KPBC_max_torque, float Joint_Bio_Sat,
                      float *Knee_torque_command, float *Ankle_torque_command,
                      float *deltaL, float *hip_pos, float *Esys, float
                      *Esys_integrate_out, float *U_S_KNEE, float *U_S_ANKLE,
                      float *COPFX, float *U_LIN_DAMP_A, float *U_STOP_K, float *
                      U_STOP_A, float *U_PBC_K, float *U_PBC_A, float
                      *knee_des_out, float *ankle_des_out, float *foot_contact,
                      float *stance, float *swing, float *phase_var_out, float
                      *IMU_LIVE_OUT, float *StanceGain, float *SwingGain, float *
                      knee_joint_vel, float *ankle_joint_vel, float *hip_vel,
                      float *PushOffOut);
extern void SLIP_KPBC_init(void);
extern void knee_pos_prev_not_empty_init(void);

#endif

/*
 * File trailer for SLIP_KPBC.h
 *
 * [EOF]
 */
