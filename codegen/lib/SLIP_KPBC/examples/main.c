/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 31-Aug-2020 12:07:22
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC.h"
#include "main.h"
#include "SLIP_KPBC_terminate.h"
#include "SLIP_KPBC_initialize.h"

/* Function Declarations */
static float argInit_real32_T(void);
static void main_SLIP_KPBC(void);

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : float
 */
static float argInit_real32_T(void)
{
  return 0.0F;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_SLIP_KPBC(void)
{
  float IMU_pitch;
  float Knee_joint_position;
  float Ankle_joint_position;
  float time_in;
  float dt;
  float Load_cell_x_force;
  float Load_cell_y_force;
  float Load_cell_z_force;
  float Load_cell_x_moment;
  float Load_cell_y_moment;
  float Load_cell_z_moment;
  float kp_knee;
  float kd_knee;
  float kp_ankle;
  float kd_ankle;
  float ankle_des_in;
  float knee_des_in;
  float vel_filter_coeff;
  float KPBC_filter_coeff;
  float SLIP_ON;
  float lt;
  float k;
  float d;
  float L0;
  float KPBC_ON;
  float KPBC_max_torque_rate;
  float pbc_gain_knee;
  float md;
  float Eref;
  float knee_stop_low;
  float Knee_torque_command;
  float Ankle_torque_command;
  float deltaL;
  float hip_pos;
  float Esys;
  float Esys_integrate_out;
  float U_S_KNEE;
  float U_S_ANKLE;
  float COPFX;
  float U_LIN_DAMP_A;
  float U_STOP_K;
  float U_STOP_A;
  float U_PBC_K;
  float U_PBC_A;
  float knee_des_out;
  float ankle_des_out;
  float foot_contact;
  float stance;
  float swing;
  float phase_var_out;
  float IMU_LIVE_OUT;
  float StanceGain;
  float SwingGain;
  float knee_joint_vel;
  float ankle_joint_vel;
  float hip_vel;
  float PushOffOut;

  /* Initialize function 'SLIP_KPBC' input arguments. */
  IMU_pitch = argInit_real32_T();
  Knee_joint_position = argInit_real32_T();
  Ankle_joint_position = argInit_real32_T();
  time_in = argInit_real32_T();
  dt = argInit_real32_T();
  Load_cell_x_force = argInit_real32_T();
  Load_cell_y_force = argInit_real32_T();
  Load_cell_z_force = argInit_real32_T();
  Load_cell_x_moment = argInit_real32_T();
  Load_cell_y_moment = argInit_real32_T();
  Load_cell_z_moment = argInit_real32_T();
  kp_knee = argInit_real32_T();
  kd_knee = argInit_real32_T();
  kp_ankle = argInit_real32_T();
  kd_ankle = argInit_real32_T();
  ankle_des_in = argInit_real32_T();
  knee_des_in = argInit_real32_T();
  vel_filter_coeff = argInit_real32_T();
  KPBC_filter_coeff = argInit_real32_T();
  SLIP_ON = argInit_real32_T();
  lt = argInit_real32_T();
  k = argInit_real32_T();
  d = argInit_real32_T();
  L0 = argInit_real32_T();
  KPBC_ON = argInit_real32_T();
  KPBC_max_torque_rate = argInit_real32_T();
  pbc_gain_knee = argInit_real32_T();
  md = argInit_real32_T();
  Eref = argInit_real32_T();
  knee_stop_low = argInit_real32_T();

  /* Call the entry-point 'SLIP_KPBC'. */
  SLIP_KPBC(IMU_pitch, Knee_joint_position, Ankle_joint_position, time_in, dt,
            Load_cell_x_force, Load_cell_y_force, Load_cell_z_force,
            Load_cell_x_moment, Load_cell_y_moment, Load_cell_z_moment, kp_knee,
            kd_knee, kp_ankle, kd_ankle, ankle_des_in, knee_des_in,
            vel_filter_coeff, KPBC_filter_coeff, SLIP_ON, lt, k, d, L0, KPBC_ON,
            KPBC_max_torque_rate, pbc_gain_knee, md, Eref, knee_stop_low,
            argInit_real32_T(), argInit_real32_T(), argInit_real32_T(),
            argInit_real32_T(), argInit_real32_T(), argInit_real32_T(),
            argInit_real32_T(), argInit_real32_T(), argInit_real32_T(),
            argInit_real32_T(), &Knee_torque_command, &Ankle_torque_command,
            &deltaL, &hip_pos, &Esys, &Esys_integrate_out, &U_S_KNEE, &U_S_ANKLE,
            &COPFX, &U_LIN_DAMP_A, &U_STOP_K, &U_STOP_A, &U_PBC_K, &U_PBC_A,
            &knee_des_out, &ankle_des_out, &foot_contact, &stance, &swing,
            &phase_var_out, &IMU_LIVE_OUT, &StanceGain, &SwingGain,
            &knee_joint_vel, &ankle_joint_vel, &hip_vel, &PushOffOut);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  SLIP_KPBC_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_SLIP_KPBC();

  /* Terminate the application.
     You do not need to do this more than one time. */
  SLIP_KPBC_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
