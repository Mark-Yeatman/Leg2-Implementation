/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: SLIP.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 15-Oct-2019 15:51:04
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP.h"
#include "sind.h"
#include "cosd.h"

/* Variable Definitions */
static float knee_pos_prev;
static boolean_T knee_pos_prev_not_empty;
static float ankle_pos_prev;
static float dknee_prev;
static float dankle_prev;
static float t;
static float t_prev;
static boolean_T t_prev_not_empty;
static float IMU_pitch_prev;

/* Function Declarations */
static float rt_powf_snf(float u0, float u1);

/* Function Definitions */

/*
 * Arguments    : float u0
 *                float u1
 * Return Type  : float
 */
static float rt_powf_snf(float u0, float u1)
{
  float y;
  float f16;
  float f17;
  if (rtIsNaNF(u0) || rtIsNaNF(u1)) {
    y = ((real32_T)rtNaN);
  } else {
    f16 = (float)fabs(u0);
    f17 = (float)fabs(u1);
    if (rtIsInfF(u1)) {
      if (f16 == 1.0F) {
        y = ((real32_T)rtNaN);
      } else if (f16 > 1.0F) {
        if (u1 > 0.0F) {
          y = ((real32_T)rtInf);
        } else {
          y = 0.0F;
        }
      } else if (u1 > 0.0F) {
        y = 0.0F;
      } else {
        y = ((real32_T)rtInf);
      }
    } else if (f17 == 0.0F) {
      y = 1.0F;
    } else if (f17 == 1.0F) {
      if (u1 > 0.0F) {
        y = u0;
      } else {
        y = 1.0F / u0;
      }
    } else if (u1 == 2.0F) {
      y = u0 * u0;
    } else if ((u1 == 0.5F) && (u0 >= 0.0F)) {
      y = (float)sqrt(u0);
    } else if ((u0 < 0.0F) && (u1 > (float)floor(u1))) {
      y = ((real32_T)rtNaN);
    } else {
      y = (float)pow(u0, u1);
    }
  }

  return y;
}

/*
 * Inputs in addition to sensors:
 * Arguments    : float IMU_pitch
 *                float Knee_motor_position
 *                float Knee_joint_position
 *                float Ankle_motor_position
 *                float Ankle_joint_position
 *                float Iteration
 *                float Iteration_time
 *                float Knee_torque_sensor
 *                float Ankle_torque_sensor
 *                float Load_cell_x_force
 *                float Load_cell_y_force
 *                float Load_cell_z_force
 *                float Load_cell_x_moment
 *                float Load_cell_y_moment
 *                float Load_cell_z_moment
 *                float kp_knee
 *                float kd_knee
 *                float kp_ankle
 *                float kd_ankle
 *                float ankle_ind_of_hip
 *                float knee_ind_of_hip
 *                float ankle_des_in
 *                float knee_des_in
 *                float time_in
 *                float filter_coeff
 *                float IMU_filter_coeff
 *                float q_h_0
 *                float q_h_min
 *                float c
 *                float s_po
 *                float FC
 *                float lf
 *                float la
 *                float ls
 *                float lt
 *                float k
 *                float d
 *                float L0
 *                float SLIP_ON
 *                float *Knee_torque_command
 *                float *Ankle_torque_command
 *                float *deltaL
 *                float *hip_pos
 *                float *t_out
 *                float *dt
 * Return Type  : void
 */
void SLIP(float IMU_pitch, float Knee_motor_position, float Knee_joint_position,
          float Ankle_motor_position, float Ankle_joint_position, float
          Iteration, float Iteration_time, float Knee_torque_sensor, float
          Ankle_torque_sensor, float Load_cell_x_force, float Load_cell_y_force,
          float Load_cell_z_force, float Load_cell_x_moment, float
          Load_cell_y_moment, float Load_cell_z_moment, float kp_knee, float
          kd_knee, float kp_ankle, float kd_ankle, float ankle_ind_of_hip, float
          knee_ind_of_hip, float ankle_des_in, float knee_des_in, float time_in,
          float filter_coeff, float IMU_filter_coeff, float q_h_0, float q_h_min,
          float c, float s_po, float FC, float lf, float la, float ls, float lt,
          float k, float d, float L0, float SLIP_ON, float *Knee_torque_command,
          float *Ankle_torque_command, float *deltaL, float *hip_pos, float
          *t_out, float *dt)
{
  boolean_T guard1 = false;
  float b_lf[4];
  boolean_T x[4];
  int b_k;
  boolean_T y;
  boolean_T exitg1;
  float F;
  float f0;
  float f1;
  float f2;
  float f3;
  float f4;
  float f5;
  float f6;
  float f7;
  float f8;
  float f9;
  float f10;
  float f11;
  float f12;
  float f13;
  float f14;
  float f15;
  float minval;
  (void)Knee_joint_position;
  (void)Ankle_joint_position;
  (void)Iteration;
  (void)Knee_torque_sensor;
  (void)Ankle_torque_sensor;
  (void)Load_cell_x_force;
  (void)Load_cell_y_force;
  (void)Load_cell_z_force;
  (void)Load_cell_x_moment;
  (void)Load_cell_y_moment;
  (void)Load_cell_z_moment;
  (void)ankle_ind_of_hip;
  (void)knee_ind_of_hip;
  (void)q_h_0;
  (void)q_h_min;
  (void)c;
  (void)s_po;
  (void)FC;

  /*  kp_ankle,kd_ankle,kp_knee,kd_knee are joint values */
  /*  ankle_ind_of_hip and knee_ind_of_hip are boolean values (true or false), */
  /*  --can be checkboxes in LabView with default TRUE. */
  /*  ankle_des_cte and knee_des_cte are desired angles for ankle and knee, */
  /*  controlled from LabView. */
  /*  Variables Defined */
  /* Persistent variables used to store data between iterations */
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
    IMU_pitch_prev = 5.0F;
  }

  /*     %% */
  /* Initialization */
  /* Knee_torque_command = 0;  */
  /* Ankle_torque_command = 0;  */
  *deltaL = 0.0F;

  /* hip_pos = 0;  */
  /* t_out = 0;  */
  /* dt = 0; */
  /*  Software/Hardware limits */
  /*      ANKLE_POS_MIN_LIM  = -35; */
  /*      ANKLE_POS_MAX_LIM = 35; */
  /*      ANKLE_VEL_MIN_LIM = -200; %deg/s */
  /*      ANKLE_VEL_MAX_LIM = 200; %deg/s */
  /*      KNEE_POS_MIN_LIM = 2; */
  /*      KNEE_POS_MAX_LIM = 85; */
  /*      KNEE_VEL_MIN_LIM = -400; %deg/s */
  /*      KNEE_VEL_MAX_LIM = 400; %deg/s */
  /* Assign joint positions */
  *hip_pos = (1.0F - IMU_filter_coeff) * IMU_pitch_prev + IMU_filter_coeff *
    IMU_pitch;

  /* Calculate joint velocities usindg weighted backwards difference */
  *dt = Iteration_time;
  t += Iteration_time;
  dknee_prev = (1.0F - filter_coeff) * dknee_prev + filter_coeff *
    (Knee_motor_position - knee_pos_prev) / Iteration_time;
  dankle_prev = (1.0F - filter_coeff) * dankle_prev + filter_coeff *
    (Ankle_motor_position - ankle_pos_prev) / Iteration_time;
  guard1 = false;
  if (SLIP_ON != 0.0F) {
    b_lf[0] = lf;
    b_lf[1] = la;
    b_lf[2] = lt;
    b_lf[3] = ls;
    for (b_k = 0; b_k < 4; b_k++) {
      x[b_k] = (b_lf[b_k] > 0.001F);
    }

    y = true;
    b_k = 0;
    exitg1 = false;
    while ((!exitg1) && (b_k < 4)) {
      if (!x[b_k]) {
        y = false;
        exitg1 = true;
      } else {
        b_k++;
      }
    }

    if (y) {
      /* Use SLIP Embedding Controller */
      /* reverse ankle and shift because sign convention of biomechanics versus biped modeling */
      F = -Ankle_motor_position;
      cosd(&F);
      f0 = Knee_motor_position;
      cosd(&f0);
      f1 = -Ankle_motor_position + Knee_motor_position;
      cosd(&f1);
      f2 = -Ankle_motor_position;
      sind(&f2);
      f3 = -Ankle_motor_position + Knee_motor_position;
      sind(&f3);
      f4 = -Ankle_motor_position;
      cosd(&f4);
      f5 = Knee_motor_position;
      cosd(&f5);
      f6 = -Ankle_motor_position + Knee_motor_position;
      cosd(&f6);
      f7 = -Ankle_motor_position;
      sind(&f7);
      f8 = -Ankle_motor_position + Knee_motor_position;
      sind(&f8);
      f9 = -Ankle_motor_position;
      cosd(&f9);
      f10 = -Ankle_motor_position + Knee_motor_position;
      cosd(&f10);
      f11 = -Ankle_motor_position;
      sind(&f11);
      f12 = -Ankle_motor_position + Knee_motor_position;
      sind(&f12);
      f13 = -Ankle_motor_position + Knee_motor_position;
      cosd(&f13);
      f14 = Knee_motor_position;
      sind(&f14);
      f15 = -Ankle_motor_position + Knee_motor_position;
      sind(&f15);
      F = -k * ((float)sqrt((((la * la + lf * lf) + ls * ls) + lt * lt) + 2.0F *
                            ((((la * ls * F + ls * lt * f0) + la * lt * f1) + lf
        * ls * f2) + lf * lt * f3)) - L0) + -d * (rt_powf_snf((((la * la + lf *
        lf) + ls * ls) + lt * lt) + 2.0F * ((((la * ls * f4 + ls * lt * f5) + la
        * lt * f6) + lf * ls * f7) + lf * lt * f8), -0.5F) * (((lf * ls * f9 +
        lf * lt * f10) + -la * (ls * f11 + lt * f12)) * -dankle_prev + lt * ((lf
        * f13 + -ls * f14) + -la * f15) * dknee_prev));
      if ((-100.0F >= F) || rtIsNaNF(F)) {
        F = -100.0F;
      }

      if ((100.0F <= F) || rtIsNaNF(F)) {
        minval = 100.0F;
      } else {
        minval = F;
      }

      F = -Ankle_motor_position;
      cosd(&F);
      f0 = -Ankle_motor_position + Knee_motor_position;
      cosd(&f0);
      f1 = -Ankle_motor_position;
      sind(&f1);
      f2 = -Ankle_motor_position + Knee_motor_position;
      sind(&f2);
      f3 = -Ankle_motor_position;
      cosd(&f3);
      f4 = Knee_motor_position;
      cosd(&f4);
      f5 = -Ankle_motor_position + Knee_motor_position;
      cosd(&f5);
      f6 = -Ankle_motor_position;
      sind(&f6);
      f7 = -Ankle_motor_position + Knee_motor_position;
      sind(&f7);
      f8 = -Ankle_motor_position + Knee_motor_position;
      cosd(&f8);
      f9 = Knee_motor_position;
      sind(&f9);
      f10 = -Ankle_motor_position + Knee_motor_position;
      sind(&f10);
      f11 = -Ankle_motor_position;
      cosd(&f11);
      f12 = Knee_motor_position;
      cosd(&f12);
      f13 = -Ankle_motor_position + Knee_motor_position;
      cosd(&f13);
      f14 = -Ankle_motor_position;
      sind(&f14);
      f15 = -Ankle_motor_position + Knee_motor_position;
      sind(&f15);
      *Ankle_torque_command = -(((lf * ls * F + lf * lt * f0) + -la * (ls * f1 +
        lt * f2)) * rt_powf_snf((((la * la + lf * lf) + ls * ls) + lt * lt) +
        2.0F * ((((la * ls * f3 + ls * lt * f4) + la * lt * f5) + lf * ls * f6)
                + lf * lt * f7), -0.5F) * minval);
      *Knee_torque_command = lt * ((lf * f8 + -ls * f9) + -la * f10) *
        rt_powf_snf((((la * la + lf * lf) + ls * ls) + lt * lt) + 2.0F * ((((la *
        ls * f11 + ls * lt * f12) + la * lt * f13) + lf * ls * f14) + lf * lt *
        f15), -0.5F) * minval;
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    /* Use PD Controller        */
    /* Calculate torque output usindddg PD controller */
    /* This control basically assumes that the desired joint angles are static i.e., the desired velocities are all zero. */
    /* function */
    /*  Helper functions */
    *Knee_torque_command = kp_knee * (knee_des_in - Knee_motor_position) +
      kd_knee * -dknee_prev;
    *Ankle_torque_command = kp_ankle * (ankle_des_in - Ankle_motor_position) +
      kd_ankle * -dankle_prev;
    if ((-150.0F >= *Knee_torque_command) || rtIsNaNF(*Knee_torque_command)) {
      F = -150.0F;
    } else {
      F = *Knee_torque_command;
    }

    if ((150.0F <= F) || rtIsNaNF(F)) {
      *Knee_torque_command = 150.0F;
    } else {
      *Knee_torque_command = F;
    }

    if ((-150.0F >= *Ankle_torque_command) || rtIsNaNF(*Ankle_torque_command)) {
      F = -150.0F;
    } else {
      F = *Ankle_torque_command;
    }

    if ((150.0F <= F) || rtIsNaNF(F)) {
      *Ankle_torque_command = 150.0F;
    } else {
      *Ankle_torque_command = F;
    }
  }

  /*     %% */
  /* Storing persistent variables for next iteration */
  knee_pos_prev = Knee_motor_position;
  ankle_pos_prev = Ankle_motor_position;
  *t_out = t;
  t_prev = time_in;
  IMU_pitch_prev = *hip_pos;
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
 * File trailer for SLIP.c
 *
 * [EOF]
 */
