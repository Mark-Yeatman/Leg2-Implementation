/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 15-Oct-2019 15:51:04
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
#include "SLIP.h"
#include "main.h"
#include "SLIP_terminate.h"
#include "SLIP_initialize.h"

/* Function Declarations */
static float argInit_real32_T(void);
static void main_SLIP(void);

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
static void main_SLIP(void)
{
  float IMU_pitch;
  float Knee_motor_position;
  float Knee_joint_position;
  float Ankle_motor_position;
  float Ankle_joint_position;
  float Iteration;
  float Iteration_time;
  float Knee_torque_sensor;
  float Ankle_torque_sensor;
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
  float ankle_ind_of_hip;
  float knee_ind_of_hip;
  float ankle_des_in;
  float knee_des_in;
  float time_in;
  float filter_coeff;
  float IMU_filter_coeff;
  float q_h_0;
  float q_h_min;
  float c;
  float Knee_torque_command;
  float Ankle_torque_command;
  float deltaL;
  float hip_pos;
  float t_out;
  float dt;

  /* Initialize function 'SLIP' input arguments. */
  IMU_pitch = argInit_real32_T();
  Knee_motor_position = argInit_real32_T();
  Knee_joint_position = argInit_real32_T();
  Ankle_motor_position = argInit_real32_T();
  Ankle_joint_position = argInit_real32_T();
  Iteration = argInit_real32_T();
  Iteration_time = argInit_real32_T();
  Knee_torque_sensor = argInit_real32_T();
  Ankle_torque_sensor = argInit_real32_T();
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
  ankle_ind_of_hip = argInit_real32_T();
  knee_ind_of_hip = argInit_real32_T();
  ankle_des_in = argInit_real32_T();
  knee_des_in = argInit_real32_T();
  time_in = argInit_real32_T();
  filter_coeff = argInit_real32_T();
  IMU_filter_coeff = argInit_real32_T();
  q_h_0 = argInit_real32_T();
  q_h_min = argInit_real32_T();
  c = argInit_real32_T();

  /* Call the entry-point 'SLIP'. */
  SLIP(IMU_pitch, Knee_motor_position, Knee_joint_position, Ankle_motor_position,
       Ankle_joint_position, Iteration, Iteration_time, Knee_torque_sensor,
       Ankle_torque_sensor, Load_cell_x_force, Load_cell_y_force,
       Load_cell_z_force, Load_cell_x_moment, Load_cell_y_moment,
       Load_cell_z_moment, kp_knee, kd_knee, kp_ankle, kd_ankle,
       ankle_ind_of_hip, knee_ind_of_hip, ankle_des_in, knee_des_in, time_in,
       filter_coeff, IMU_filter_coeff, q_h_0, q_h_min, c, argInit_real32_T(),
       argInit_real32_T(), argInit_real32_T(), argInit_real32_T(),
       argInit_real32_T(), argInit_real32_T(), argInit_real32_T(),
       argInit_real32_T(), argInit_real32_T(), argInit_real32_T(),
       &Knee_torque_command, &Ankle_torque_command, &deltaL, &hip_pos, &t_out,
       &dt);
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
  SLIP_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_SLIP();

  /* Terminate the application.
     You do not need to do this more than one time. */
  SLIP_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
