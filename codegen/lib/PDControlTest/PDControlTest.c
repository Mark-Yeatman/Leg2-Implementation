/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: PDControlTest.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 03-Oct-2019 15:42:50
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PDControlTest.h"

/* Variable Definitions */
static float knee_pos_prev;
static boolean_T knee_pos_prev_not_empty;
static float ankle_pos_prev;
static float dknee_prev;
static float dankle_prev;
static float t;
static float t_prev;
static boolean_T t_prev_not_empty;

/* Function Definitions */

/*
 * Variables Defined
 * Persistent variables used to store data between iterations
 * Arguments    : float knee_pos
 *                float ankle_pos
 *                float knee_des
 *                float ankle_des
 *                float dt
 *                float kp_k
 *                float kd_k
 *                float kp_a
 *                float kd_a
 *                float time_in
 *                float filter_coeff
 *                float *u_k
 *                float *u_a
 *                float *t_out
 *                float *dt_out
 * Return Type  : void
 */
void PDControlTest(float knee_pos, float ankle_pos, float knee_des, float
                   ankle_des, float dt, float kp_k, float kd_k, float kp_a,
                   float kd_a, float time_in, float filter_coeff, float *u_k,
                   float *u_a, float *t_out, float *dt_out)
{
  float maxval;
  if (!t_prev_not_empty) {
    t_prev = time_in;
    t_prev_not_empty = true;
  }

  if ((!knee_pos_prev_not_empty) || (time_in - t_prev > 1.0F)) {
    knee_pos_prev = 0.0F;
    knee_pos_prev_not_empty = true;
    ankle_pos_prev = 0.0F;
    t = 0.0F;
    dknee_prev = 0.0F;
    dankle_prev = 0.0F;
  }

  /*  Initialization */
  /*  Software/Hardware limits */
  /*      ANKLE_POS_MIN_LIM  = -35; */
  /*      ANKLE_POS_MAX_LIM = 35; */
  /*      ANKLE_VEL_MIN_LIM = -200; %deg/s */
  /*      ANKLE_VEL_MAX_LIM = 200; %deg/s */
  /*      KNEE_POS_MIN_LIM = 2; */
  /*      KNEE_POS_MAX_LIM = 85; */
  /*      KNEE_VEL_MIN_LIM = -400; %deg/s */
  /*      KNEE_VEL_MAX_LIM = 400; %deg/s */
  /*  Control */
  /* Calculate joint velocities using weighted backwards difference */
  t += dt;
  dknee_prev = (1.0F - filter_coeff) * dknee_prev + filter_coeff * (knee_pos -
    knee_pos_prev) / dt;
  dankle_prev = (1.0F - filter_coeff) * dankle_prev + filter_coeff * (ankle_pos
    - ankle_pos_prev) / dt;

  /* [u_k,u_a] = PDControl(kp_k,kd_k,knee_des,knee_pos,dknee,kp_a,kd_a,ankle_des,ankle_pos,dankle); */
  *u_k = kp_k * (knee_des - knee_pos) + kd_k * -dknee_prev;
  *u_a = kp_a * (ankle_des - ankle_pos) + kd_a * -dankle_prev;
  if ((-60.0F >= *u_k) || rtIsNaNF(*u_k)) {
    maxval = -60.0F;
  } else {
    maxval = *u_k;
  }

  if ((60.0F <= maxval) || rtIsNaNF(maxval)) {
    *u_k = 60.0F;
  } else {
    *u_k = maxval;
  }

  if ((-60.0F >= *u_a) || rtIsNaNF(*u_a)) {
    maxval = -60.0F;
  } else {
    maxval = *u_a;
  }

  if ((60.0F <= maxval) || rtIsNaNF(maxval)) {
    *u_a = 60.0F;
  } else {
    *u_a = maxval;
  }

  /*  Storing persistent variables for next iteration */
  knee_pos_prev = knee_pos;
  ankle_pos_prev = ankle_pos;
  *t_out = t;
  t_prev = time_in;
  *dt_out = dt;

  /*  %% Helper functions */
  /*  function [Knee_torque_command,Ankle_torque_command] = PDControl(kp_k, kd_k, q_kstar, q_k, qdot_k, kp_a, kd_a, q_astar, q_a, qdot_a) */
  /*      %This control basically assumes that the desired joint angles are static i.e., the desired velocities are all zero. */
  /*      Knee_torque_command = kp_k*(q_kstar-q_k) + kd_k*(-qdot_k); */
  /*      Ankle_torque_command = kp_a*(q_astar-q_a) + kd_a*(-qdot_a); */
  /*  end */
  /*   */
  /*  function y = Saturate(x,x1,x2) */
  /*      %Function to prevent the desired joint angles from changing to fast.  */
  /*      %Works via saturation */
  /*      y=min(x,max(x1,x2)); */
  /*      y=max(y,min(x1,x2)); */
  /*  end */
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void knee_pos_prev_not_empty_init(void)
{
  knee_pos_prev_not_empty = false;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void t_prev_not_empty_init(void)
{
  t_prev_not_empty = false;
}

/*
 * File trailer for PDControlTest.c
 *
 * [EOF]
 */
