/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: SLIP_KPBC.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 01-Jul-2020 11:22:48
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC.h"
#include "sind.h"
#include "cosd.h"
#include "SLIP_KPBC_rtwutil.h"

/* Variable Definitions */
static float hip_pos_prev;
static float knee_pos_prev;
static boolean_T knee_pos_prev_not_empty;
static float ankle_pos_prev;
static float knee_vel_prev;
static float ankle_vel_prev;
static float hip_vel_prev;
static float t_switched_phase;
static boolean_T Stance;
static boolean_T Swing;
static float ForceCount;
static float u_pbc_knee_prev;
static float u_pbc_ankle_prev;

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
  float f1;
  float f2;
  if (rtIsNaNF(u0) || rtIsNaNF(u1)) {
    y = ((real32_T)rtNaN);
  } else {
    f1 = (float)fabs(u0);
    f2 = (float)fabs(u1);
    if (rtIsInfF(u1)) {
      if (f1 == 1.0F) {
        y = ((real32_T)rtNaN);
      } else if (f1 > 1.0F) {
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
    } else if (f2 == 0.0F) {
      y = 1.0F;
    } else if (f2 == 1.0F) {
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
 *                float Knee_joint_position
 *                float Ankle_joint_position
 *                float time_in
 *                float dt
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
 *                float ankle_des_in
 *                float knee_des_in
 *                float vel_filter_coeff
 *                float KPBC_filter_coeff
 *                float SLIP_ON
 *                float lt
 *                float k
 *                float d
 *                float L0
 *                float KPBC_ON
 *                float KPBC_max_torque_rate
 *                float pbc_gain_knee
 *                float M
 *                float Eref
 *                float knee_stop_low
 *                float knee_stop_high
 *                float ankle_stop_low
 *                float ankle_stop_high
 *                float max_torque
 *                float k_tilde
 *                float Ignore_PushOff
 *                float F_thresh
 *                float Command_State
 *                float KPBC_max_torque
 *                float Joint_Bio_Sat
 *                float *Knee_torque_command
 *                float *Ankle_torque_command
 *                float *deltaL
 *                float *hip_pos
 *                float *Esys
 *                float *Esys_integrate_out
 *                float *U_LIN_SPRING_K
 *                float *U_LIN_SPRING_A
 *                float *U_LIN_DAMP_K
 *                float *U_LIN_DAMP_A
 *                float *U_STOP_K
 *                float *U_STOP_A
 *                float *U_PBC_K
 *                float *U_PBC_A
 *                float *knee_des_out
 *                float *ankle_des_out
 *                float *foot_contact
 *                float *stance
 *                float *swing
 *                float *phase_var_out
 *                float *IMU_LIVE_OUT
 *                float *StanceGain
 *                float *SwingGain
 *                float *knee_joint_vel
 *                float *ankle_joint_vel
 *                float *hip_vel
 *                float *PushOffOut
 * Return Type  : void
 */
void SLIP_KPBC(float IMU_pitch, float Knee_joint_position, float
               Ankle_joint_position, float time_in, float dt, float
               Load_cell_x_force, float Load_cell_y_force, float
               Load_cell_z_force, float Load_cell_x_moment, float
               Load_cell_y_moment, float Load_cell_z_moment, float kp_knee,
               float kd_knee, float kp_ankle, float kd_ankle, float ankle_des_in,
               float knee_des_in, float vel_filter_coeff, float
               KPBC_filter_coeff, float SLIP_ON, float lt, float k, float d,
               float L0, float KPBC_ON, float KPBC_max_torque_rate, float
               pbc_gain_knee, float M, float Eref, float knee_stop_low, float
               knee_stop_high, float ankle_stop_low, float ankle_stop_high,
               float max_torque, float k_tilde, float Ignore_PushOff, float
               F_thresh, float Command_State, float KPBC_max_torque, float
               Joint_Bio_Sat, float *Knee_torque_command, float
               *Ankle_torque_command, float *deltaL, float *hip_pos, float *Esys,
               float *Esys_integrate_out, float *U_LIN_SPRING_K, float
               *U_LIN_SPRING_A, float *U_LIN_DAMP_K, float *U_LIN_DAMP_A, float *
               U_STOP_K, float *U_STOP_A, float *U_PBC_K, float *U_PBC_A, float *
               knee_des_out, float *ankle_des_out, float *foot_contact, float
               *stance, float *swing, float *phase_var_out, float *IMU_LIVE_OUT,
               float *StanceGain, float *SwingGain, float *knee_joint_vel, float
               *ankle_joint_vel, float *hip_vel, float *PushOffOut)
{
  float t7;
  float f0;
  float t8;
  float t11;
  float t14;
  float value;
  float t19;
  float t20;
  float t23;
  float t27;
  float t24;
  float t25;
  float t26;
  float t28;
  float t33;
  float t29;
  float t30;
  float J[5];
  boolean_T IMU_LIVE;
  float x[3];
  int i;
  boolean_T FootContact;
  signed char i0;
  float u_lin_spring[5];
  float u_lin_damp[5];
  float u_stop[2];
  float knee_limits[2];
  float ankle_limits[2];
  int trueCount;
  static const short iv0[8] = { 1, 139, 400, 530, 720, 848, 976, 1001 };

  static const short knee_ind[8] = { 1, 139, 400, 530, 720, 848, 976, 1001 };

  static const float a_knee[28] = { -21.9472504F, 20.9295F, 15.9362497F,
    -59.8412514F, -19.4405F, 45.107F, 1.42525F, 4.24525F, -6.9765F, -0.99825F,
    13.6422501F, -1.4025F, -7.153F, 0.06575F, 21.672F, 7.719F, 18.664F, 64.863F,
    38.41F, 0.456F, 2.21F, 0.0F, 0.0F, 25.883F, 0.0F, -47.296F, 0.0F, 3.245F };

  static const short iv1[11] = { 1, 61, 205, 330, 440, 550, 653, 846, 898, 977,
    1001 };

  static const short ankle_ind[11] = { 1, 61, 205, 330, 440, 550, 653, 846, 898,
    977, 1001 };

  static const float a_ankle[40] = { 5.5005F, -8.985F, 0.8125F, -2.19975F,
    -14.7575F, 23.0932F, -29.802F, 0.7035F, -2.6085F, -0.729F, -0.8805F, 3.7702F,
    -0.475F, 0.32075F, 0.8785F, -3.9142F, 9.934F, -0.2345F, 0.8695F, 0.019F,
    -4.6F, 5.266F, 7.741F, 9.62F, -0.745F, -19.924F, -0.056F, -0.525F, 1.214F,
    0.58F, 0.0F, 4.6512F, 2.8125F, 0.0F, -24.244F, 0.0F, 0.0F, 0.0F, 0.0F,
    -1.344F };

  float b_t20;
  float b_t7;
  float b_t26;
  float c_t7;
  (void)Load_cell_x_moment;
  (void)Load_cell_y_moment;
  (void)Load_cell_z_moment;
  (void)Ignore_PushOff;

  /*  Knee_joint_position, Ankle_joint_position are in degrees, in a */
  /*  "biomechanics frame" */
  /*  t, dt, t_out are in seconds */
  /*     %% Persistent Variable Definitions */
  /* Persistent variables used to store data between iterations */
  /* t_switched_swing... %implement swing setpoint switching hold? */
  if (!knee_pos_prev_not_empty) {
    hip_pos_prev = IMU_pitch;
    knee_pos_prev_not_empty = true;
    t_switched_phase = time_in;

    /* t_switched_swing = time_in; */
  }

  /*     %% Initialization and Hard Coded Values */
  /* Knee_torque_command = 0;  */
  /* Ankle_torque_command = 0;  */
  /* deltaL = single(0); %meters */
  /* Esys = single(0); %Joules */
  /* g = single(9.81); %m/s^2 */
  /* pulse_length = single(10); %seconds */
  /* hip_pos = 0;  */
  /* t_out = 0;  */
  /* dt = 0; */
  /*  converted to kg */
  /* carbon fiber foot + mechanism, kg */
  /* lt = single(0.3733); %meters */
  /* meters */
  /* meters */
  /* meters */
  /* meters */
  /* meters */
  /* From ordering in makeMatlabFunctionsProthesisTestBench */
  *U_LIN_DAMP_K = 0.0F;
  *U_LIN_DAMP_A = 0.0F;
  *U_LIN_SPRING_K = 0.0F;
  *U_LIN_SPRING_A = 0.0F;
  *U_PBC_K = 0.0F;
  *U_PBC_A = 0.0F;

  /*  Position/velocity hard limits */
  /* deg/s */
  /* deg/s */
  /* deg/s */
  /* deg/s */
  /* joules */
  /* Winters Data Arrays */
  /*     %% State Calculations */
  /* Assign joint positions */
  *hip_pos = IMU_pitch;

  /* Calculate joint velocities usindg exponential smoothing filter */
  /* https://en.wikipedia.org/wiki/Exponential_smoothing */
  *hip_vel = (1.0F - vel_filter_coeff) * hip_vel_prev + vel_filter_coeff *
    (IMU_pitch - hip_pos_prev) / dt;
  knee_vel_prev = (1.0F - vel_filter_coeff) * knee_vel_prev + vel_filter_coeff *
    (Knee_joint_position - knee_pos_prev) / dt;
  ankle_vel_prev = (1.0F - vel_filter_coeff) * ankle_vel_prev + vel_filter_coeff
    * (Ankle_joint_position - ankle_pos_prev) / dt;

  /* Human Leg as a Robot states */
  /* x,y,-knee ankle, ankle angle */
  /* reverse knee sign convention of biomechanics versus biped modeling      */
  /* reverse knee sign convention of biomechanics versus biped modeling */
  /* Virtual linear spring from hip to foot, l1 is shank length, l2 is thigh length      */
  /* SPRING_JACOBIAN_FUNC */
  /*     LJACOB = SPRING_JACOBIAN_FUNC(IN1,IN2) */
  /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
  /*     31-Jan-2020 14:49:45 */
  /* Prosthesis file. Needs state and parameters as inputs */
  t7 = -Knee_joint_position + Ankle_joint_position;
  f0 = 0.0F;
  cosd(&f0);
  t8 = lt * f0;
  f0 = 0.0F;
  sind(&f0);
  t11 = lt * f0;
  t14 = -Knee_joint_position + Ankle_joint_position;
  value = t14;
  cosd(&value);
  sind(&t14);
  f0 = -Knee_joint_position;
  cosd(&f0);
  t19 = 0.3733F * f0;
  f0 = -Knee_joint_position;
  sind(&f0);
  t20 = 0.3733F * f0;
  f0 = t7;
  cosd(&f0);
  t23 = 0.15F * lt * f0;
  sind(&t7);
  t27 = 0.0628F * lt * t7;
  t24 = 0.0F * t19;
  t25 = 0.0628F * value;
  t26 = 0.15F * value;
  t28 = 0.0F * t20;
  t33 = t14 * 0.0F;
  t29 = 0.0F * t25;
  t30 = 0.0F * t26;
  t7 = ((t11 + t20) + 0.0628F * t14) + t26;
  t20 = ((t8 + t19) + t25) + -(0.15F * t14);
  t14 = 1.0F / (float)sqrt(t20 * t20 + t7 * t7);
  f0 = -Knee_joint_position;
  sind(&f0);
  t20 = Ankle_joint_position;
  cosd(&t20);
  t26 = Ankle_joint_position;
  sind(&t26);
  J[0] = 0.0F;
  J[1] = 0.0F;
  J[2] = -t14 * ((((((t24 + t28) + t29) + t30) + t33) + 0.0F * t8) + 0.0F * t11);
  J[3] = -t14 * ((((((t23 + t24) + t27) + t28) + t33) + value * 0.0F) + 0.3733F *
                 lt * f0);
  J[4] = -t14 * ((((((t23 + t27) + t29) + t30) + t33) + 0.055995F * t20) +
                 0.0234432388F * t26);

  /* SPRING_LENGTH_FUNC */
  /*     L1 = SPRING_LENGTH_FUNC(IN1,IN2) */
  /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
  /*     31-Jan-2020 14:49:42 */
  /* Prosthesis file. Needs state and parameters as inputs */
  t7 = -Knee_joint_position + Ankle_joint_position;
  value = t7;
  cosd(&value);
  sind(&t7);
  f0 = -Knee_joint_position;
  sind(&f0);
  t20 = 0.0F;
  sind(&t20);
  t24 = ((0.0628F * t7 + 0.15F * value) + 0.3733F * f0) + lt * t20;
  f0 = -Knee_joint_position;
  cosd(&f0);
  t20 = 0.0F;
  cosd(&t20);
  t19 = ((0.0628F * value - 0.15F * t7) + 0.3733F * f0) + lt * t20;
  *deltaL = (float)sqrt(t24 * t24 + t19 * t19) - L0;

  /* SPRING_VEL_FUNC */
  /*     L1DOT = SPRING_VEL_FUNC(IN1,IN2) */
  /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
  /*     31-Jan-2020 14:49:44 */
  /* Prosthesis file. Needs state and parameters as inputs */
  value = 0.0F;
  cosd(&value);
  t20 = 0.0F;
  sind(&t20);
  t14 = -Knee_joint_position;
  cosd(&t14);
  t26 = -Knee_joint_position;
  sind(&t26);
  t8 = -Knee_joint_position + Ankle_joint_position;
  t7 = -knee_vel_prev + ankle_vel_prev;
  t25 = t8;
  cosd(&t25);
  sind(&t8);
  t24 = ((0.0628F * t25 - 0.15F * t8) + 0.3733F * t14) + lt * value;
  t19 = ((0.0628F * t8 + 0.15F * t25) + 0.3733F * t26) + lt * t20;
  t25 = 1.0F / (float)sqrt(t24 * t24 + t19 * t19) * ((((0.0628F * t7 * t25 -
    0.15F * t7 * t8) + 0.3733F * -knee_vel_prev * t14) + lt * value * 0.0F) *
    (((0.0628F * t8 * 2.0F + 0.15F * t25 * 2.0F) + 0.3733F * t26 * 2.0F) + lt *
     t20 * 2.0F) - (((0.0628F * t7 * t8 + 0.15F * t7 * t25) + 0.3733F *
                     -knee_vel_prev * t26) + lt * t20 * 0.0F) * (((0.0628F * t25
    * 2.0F - 0.15F * t8 * 2.0F) + 0.3733F * t14 * 2.0F) + lt * value * 2.0F)) /
    2.0F;

  /* Calculate system energy */
  /* Ehip = 1/2*M*(Ldot^2) + (L-L0)*M*g; %kinetic and potential energy of slip model */
  /* virtual spring potential energy */
  /* KE_func(x,params); */
  /* PE = PE_func(x,params); */
  *Esys = 0.5F * k_tilde * (*deltaL * *deltaL) + 0.5F * M * (t25 * t25);

  /* Espring + PE + KE; */
  /* Phase Variable  */
  /* s_a=clamp(1+(1-s_m)/(q_h_0-q_h_m)*(hip-q_h_0),0,1); */
  t24 = 0.4F * (IMU_pitch - -20.0F) / 43.0F;
  if (0.6F + t24 <= 1.0F) {
    t26 = 0.6F + t24;
  } else {
    t26 = 1.0F;
  }

  if (t26 >= 0.6F) {
    *phase_var_out = t26;
  } else {
    *phase_var_out = 0.6F;
  }

  /* IMU State check */
  IMU_LIVE = ((float)fabs(hip_pos_prev - IMU_pitch) < 1.0F);

  /* PushOff = and(or((hip_pos < hip_thresh_st),PushOff),Stance); %PushOff is persistent, so this remembers when the hip_pos crosses and stays high until Stance is off */
  /*     %% Control logic */
  /* Determine foot contact */
  x[0] = Load_cell_x_force;
  x[1] = Load_cell_y_force;
  x[2] = Load_cell_z_force;
  t7 = 0.0F;
  t20 = 1.17549435E-38F;
  for (i = 0; i < 3; i++) {
    t14 = (float)fabs(x[i]);
    if (t14 > t20) {
      t26 = t20 / t14;
      t7 = 1.0F + t7 * t26 * t26;
      t20 = t14;
    } else {
      t26 = t14 / t20;
      t7 += t26 * t26;
    }
  }

  t7 = t20 * (float)sqrt(t7);
  if (t7 > F_thresh) {
    ForceCount++;
    FootContact = (ForceCount > 2.0F);
  } else {
    FootContact = false;
    ForceCount = 0.0F;
  }

  *StanceGain = 0.0F;
  *SwingGain = 0.0F;

  /* Phase State management, set autodetection or command stance/swing */
  f0 = rt_roundf_snf(Command_State);
  if (f0 < 128.0F) {
    if (f0 >= -128.0F) {
      i0 = (signed char)f0;
    } else {
      i0 = MIN_int8_T;
    }
  } else if (f0 >= 128.0F) {
    i0 = MAX_int8_T;
  } else {
    i0 = 0;
  }

  if (i0 == 0) {
    /* Determine stance vs swing, only allow phase state switching every t_hold seconds */
    if (Stance && (!FootContact)) {
      /* The state needs to be switched */
      if (0.3F < time_in - t_switched_phase) {
        /* Check the last time we switched, if its too fast, wait */
        /* Switch to SWING */
        t_switched_phase = time_in;
        Swing = true;
        Stance = false;

        /* reset the energy intergration variables */
        /* Esys_integrate = 0; */
        /* Edis = 0;  */
        u_pbc_knee_prev = 0.0F;
        u_pbc_ankle_prev = 0.0F;
      }
    } else {
      if (Swing && FootContact && (0.3F < time_in - t_switched_phase)) {
        /* ditto, but for the other case */
        /* Switch to STANCE */
        t_switched_phase = time_in;
        Swing = false;
        Stance = true;
      }
    }

    /* Manage torque smoothing when switching phases */
    if (Stance) {
      t7 = (time_in - t_switched_phase) / 0.24000001F;
      if (t7 <= 1.0F) {
        *StanceGain = t7;
      } else {
        *StanceGain = 1.0F;
      }

      *SwingGain = 1.0F - *StanceGain;
    } else {
      if (Swing) {
        t7 = (time_in - t_switched_phase) / 0.24000001F;
        if (t7 <= 1.0F) {
          *SwingGain = t7;
        } else {
          *SwingGain = 1.0F;
        }

        *StanceGain = 1.0F - *SwingGain;
      }
    }
  } else {
    f0 = rt_roundf_snf(Command_State);
    if (f0 < 128.0F) {
      if (f0 >= -128.0F) {
        i0 = (signed char)f0;
      } else {
        i0 = MIN_int8_T;
      }
    } else if (f0 >= 128.0F) {
      i0 = MAX_int8_T;
    } else {
      i0 = 0;
    }

    if (i0 == 1) {
      /* Stance */
      Stance = true;
      Swing = false;
      *StanceGain = 1.0F;
    } else {
      f0 = rt_roundf_snf(Command_State);
      if (f0 < 128.0F) {
        if (f0 >= -128.0F) {
          i0 = (signed char)f0;
        } else {
          i0 = MIN_int8_T;
        }
      } else if (f0 >= 128.0F) {
        i0 = MAX_int8_T;
      } else {
        i0 = 0;
      }

      if (i0 == 2) {
        /* Swing */
        Stance = false;
        Swing = true;
        *SwingGain = 1.0F;
      }
    }
  }

  /* Torque computation */
  if (SLIP_ON != 0.0F) {
    /*         %% Stance                 */
    /* Virtual linear spring */
    t24 = -k * *deltaL;
    t19 = -d * t25;
    for (i = 0; i < 5; i++) {
      u_lin_spring[i] = t24 * J[i];
      u_lin_damp[i] = t19 * J[i];
    }

    *U_LIN_DAMP_K = -u_lin_damp[3];

    /*  -u_G_shape(4); */
    *U_LIN_DAMP_A = u_lin_damp[4];

    /*  u_G_shape(5); */
    *U_LIN_SPRING_K = -u_lin_spring[3];
    *U_LIN_SPRING_A = u_lin_spring[4];
    t19 = -u_lin_damp[3] + -u_lin_spring[3];

    /* reverse knee sign convention of biomechanics versus biped modeling */
    t27 = u_lin_damp[4] + u_lin_spring[4];

    /* estimate how much energy the virtual linear damper is dissipating, to be added back in during pushoff */
    /*          if ~PushOff */
    /*              Edis = Edis  + (knee_vel*u_lin_damp(4) + ankle_vel*u_lin_damp(5))*dt; */
    /*          end */
    if (KPBC_ON != 0.0F) {
      /* Also use energy tracking controller    */
      /*              if Joint_Bio_Sat  */
      /*                  %Joint level KPBC */
      /*                  u_pbc(4) = -pbc_gain_knee * knee_vel * (Esys_integrate - Eref - max(Edis,Edismax)); */
      /*                  u_pbc(5) = -pbc_gain_ankle * ankle_vel * (Esys_integrate - Eref - max(Edis,Edismax)); */
      /*                   */
      /*                  U_PBC_K = u_pbc(4); %THE SIGN NEEDS TO BE CHANGED IF USING SPRING VERSUS JOINT LEVEL BECAUSE OF BIO VS ROBOT KNEE CORDS */
      /*                  U_PBC_A = u_pbc(5); */
      /*                   */
      /*              else */
      /* KPBC along spring axis */
      t24 = *Esys - Eref;
      for (i = 0; i < 5; i++) {
        J[i] = -J[i] * pbc_gain_knee * t24 * t25;
      }

      *U_PBC_K = -J[3];

      /* THE SIGN NEEDS TO BE CHANGED IF USING SPRING VERSUS JOINT LEVEL BECAUSE OF BIO VS ROBOT KNEE CORDS */
      *U_PBC_A = J[4];

      /* end */
      /* moving average filter and saturation laws on value and rate of */
      /* torque */
      /* Absolute saturation */
      if (KPBC_max_torque > 0.0F) {
        /* Function to prevent the desired joint angles from changing to fast.  */
        /* Works via saturation */
        if (-KPBC_max_torque >= KPBC_max_torque) {
          t14 = -KPBC_max_torque;
        } else {
          t14 = KPBC_max_torque;
        }

        if (-J[3] <= t14) {
          t26 = -J[3];
        } else {
          t26 = t14;
        }

        if (t26 >= -KPBC_max_torque) {
          *U_PBC_K = t26;
        } else {
          *U_PBC_K = -KPBC_max_torque;
        }

        /* Function to prevent the desired joint angles from changing to fast.  */
        /* Works via saturation */
        if (-KPBC_max_torque >= KPBC_max_torque) {
          t14 = -KPBC_max_torque;
        } else {
          t14 = KPBC_max_torque;
        }

        if (J[4] <= t14) {
          t26 = J[4];
        } else {
          t26 = t14;
        }

        if (t26 >= -KPBC_max_torque) {
          *U_PBC_A = t26;
        } else {
          *U_PBC_A = -KPBC_max_torque;
        }
      }

      /* Rate saturation */
      if (KPBC_max_torque_rate > 0.0F) {
        t7 = (*U_PBC_K - u_pbc_knee_prev) / dt;
        t20 = (*U_PBC_A - u_pbc_ankle_prev) / dt;

        /* Function to prevent the desired joint angles from changing to fast.  */
        /* Works via saturation */
        if (-KPBC_max_torque_rate >= KPBC_max_torque_rate) {
          t14 = -KPBC_max_torque_rate;
        } else {
          t14 = KPBC_max_torque_rate;
        }

        if (t7 <= t14) {
          t26 = t7;
        } else {
          t26 = t14;
        }

        /* Function to prevent the desired joint angles from changing to fast.  */
        /* Works via saturation */
        if (-KPBC_max_torque_rate >= KPBC_max_torque_rate) {
          t14 = -KPBC_max_torque_rate;
        } else {
          t14 = KPBC_max_torque_rate;
        }

        if (t20 <= t14) {
          t7 = t20;
        } else {
          t7 = t14;
        }

        if (t26 >= -KPBC_max_torque_rate) {
          b_t26 = t26;
        } else {
          b_t26 = -KPBC_max_torque_rate;
        }

        *U_PBC_K = b_t26 * dt + u_pbc_knee_prev;
        if (t7 >= -KPBC_max_torque_rate) {
          c_t7 = t7;
        } else {
          c_t7 = -KPBC_max_torque_rate;
        }

        *U_PBC_A = c_t7 * dt + u_pbc_ankle_prev;
      }

      /* exponential smoothing of torque */
      *U_PBC_K = (1.0F - KPBC_filter_coeff) * u_pbc_knee_prev +
        KPBC_filter_coeff * *U_PBC_K;
      *U_PBC_A = (1.0F - KPBC_filter_coeff) * u_pbc_ankle_prev +
        KPBC_filter_coeff * *U_PBC_A;
      if (Joint_Bio_Sat != 0.0F) {
        /* Biomemetic power saturation */
        t7 = knee_vel_prev * *U_PBC_K;
        t20 = ankle_vel_prev * *U_PBC_A;
        if (t20 < 0.0F) {
          b_t20 = -1.0F;
        } else if (t20 > 0.0F) {
          b_t20 = 1.0F;
        } else if (t20 == 0.0F) {
          b_t20 = 0.0F;
        } else {
          b_t20 = t20;
        }

        if (b_t20 == -1.0F) {
          *U_PBC_A = 0.0F;
        }

        if (t7 < 0.0F) {
          b_t7 = -1.0F;
        } else if (t7 > 0.0F) {
          b_t7 = 1.0F;
        } else if (t7 == 0.0F) {
          b_t7 = 0.0F;
        } else {
          b_t7 = t7;
        }

        if (b_t7 == 1.0F) {
          *U_PBC_K = 0.0F;
        }
      }

      /* Use energy intergration scheme   */
      /* Esys_integrate = Esys_integrate  + (knee_vel*U_PBC_K + ankle_vel*U_PBC_A)*dt; */
      /* cancel damping term during KPBC */
      t19 = (t19 + *U_PBC_K) - (-u_lin_damp[3]);
      t27 = (t27 + *U_PBC_A) - u_lin_damp[4];
    }

    /*         %% Swing */
    if (IMU_LIVE) {
      /* Use hip pos as phase variable to index winters data */
      t24 = *phase_var_out * 1000.0F;
      trueCount = 0;
      for (i = 0; i < 8; i++) {
        if (t24 - (float)knee_ind[i] >= 0.0F) {
          trueCount++;
        }
      }

      t7 = (t24 - (float)iv0[(int)(float)trueCount - 1]) / (float)(iv0[(int)
        ((float)trueCount + 1.0F) - 1] - iv0[(int)(float)trueCount - 1]);
      *knee_des_out = ((a_knee[(int)(float)trueCount - 1] * ((t7 - 1.0F) * (t7 -
        1.0F)) + a_knee[(int)(float)trueCount + 6] * rt_powf_snf(t7 - 1.0F, 6.0F))
                       + a_knee[(int)(float)trueCount + 20] * (t7 - 1.0F)) +
        a_knee[(int)(float)trueCount + 13];
      trueCount = 0;
      for (i = 0; i < 11; i++) {
        if (t24 - (float)ankle_ind[i] >= 0.0F) {
          trueCount++;
        }
      }

      t7 = (t24 - (float)iv1[(int)(float)trueCount - 1]) / (float)(iv1[(int)
        ((float)trueCount + 1.0F) - 1] - iv1[(int)(float)trueCount - 1]);
      *ankle_des_out = ((a_ankle[(int)(float)trueCount - 1] * ((t7 - 1.0F) * (t7
        - 1.0F)) + a_ankle[(int)(float)trueCount + 9] * rt_powf_snf(t7 - 1.0F,
        6.0F)) + a_ankle[(int)(float)trueCount + 29] * (t7 - 1.0F)) + a_ankle
        [(int)(float)trueCount + 19];

      /* follow the trjactory via PD controller */
      /* This control basically assumes that the desired joint angles are static i.e., the desired velocities are all zero. */
      /*  Helper functions */
      t7 = kp_knee * (*knee_des_out - Knee_joint_position) + kd_knee *
        -knee_vel_prev;
      t20 = kp_ankle * (*ankle_des_out - Ankle_joint_position) + kd_ankle *
        -ankle_vel_prev;
    } else {
      *knee_des_out = Knee_joint_position;
      *ankle_des_out = Ankle_joint_position;

      /* This control basically assumes that the desired joint angles are static i.e., the desired velocities are all zero. */
      /*  Helper functions */
      t7 = kp_knee * (Knee_joint_position - Knee_joint_position) + kd_knee *
        -knee_vel_prev;
      t20 = kp_ankle * (Ankle_joint_position - Ankle_joint_position) + kd_ankle *
        -ankle_vel_prev;
    }

    /* Smooth between stance and swing torques */
    *Knee_torque_command = *StanceGain * t19 + *SwingGain * t7;
    *Ankle_torque_command = *StanceGain * t27 + *SwingGain * t20;
  } else {
    IMU_LIVE = true;

    /* Use Setpoint PD control */
    *knee_des_out = knee_des_in;
    *ankle_des_out = ankle_des_in;

    /* This control basically assumes that the desired joint angles are static i.e., the desired velocities are all zero. */
    /*  Helper functions */
    *Knee_torque_command = kp_knee * (knee_des_in - Knee_joint_position) +
      kd_knee * -knee_vel_prev;
    *Ankle_torque_command = kp_ankle * (ankle_des_in - Ankle_joint_position) +
      kd_ankle * -ankle_vel_prev;
  }

  /*     %% Virtual hard stops */
  for (i = 0; i < 2; i++) {
    u_stop[i] = 0.0F;
  }

  if ((knee_stop_low < knee_stop_high) && (ankle_stop_low < ankle_stop_high)) {
    knee_limits[0] = knee_stop_low;
    knee_limits[1] = knee_stop_high;

    /* deg, need to make sure its ordered */
    ankle_limits[0] = ankle_stop_low;
    ankle_limits[1] = ankle_stop_high;
  } else {
    /* deg, need to make sure its ordered */
    for (i = 0; i < 2; i++) {
      knee_limits[i] = 2.0F + 103.0F * (float)i;
      ankle_limits[i] = -35.0F + 70.0F * (float)i;
    }
  }

  if (Knee_joint_position < knee_limits[0]) {
    u_stop[0] = -kp_knee * (Knee_joint_position - knee_limits[0]) - kd_knee *
      knee_vel_prev;
  }

  if (knee_limits[1] < Knee_joint_position) {
    u_stop[0] = -kp_knee * (Knee_joint_position - knee_limits[1]) - kd_knee *
      knee_vel_prev;
  }

  if (Ankle_joint_position < ankle_limits[0]) {
    u_stop[1] = -kp_ankle * (Ankle_joint_position - ankle_limits[0]) - kd_ankle *
      ankle_vel_prev;
  }

  if (ankle_limits[1] < Ankle_joint_position) {
    u_stop[1] = -kp_ankle * (Ankle_joint_position - ankle_limits[1]) - kd_ankle *
      ankle_vel_prev;
  }

  *U_STOP_K = u_stop[0];
  *U_STOP_A = u_stop[1];
  *Knee_torque_command += u_stop[0];
  *Ankle_torque_command += u_stop[1];

  /*     %% Saturate torque output */
  if (max_torque > 0.0F) {
    /* Function to prevent the desired joint angles from changing to fast.  */
    /* Works via saturation */
    if (-max_torque >= max_torque) {
      t14 = -max_torque;
    } else {
      t14 = max_torque;
    }

    if (*Knee_torque_command <= t14) {
      t26 = *Knee_torque_command;
    } else {
      t26 = t14;
    }

    if (t26 >= -max_torque) {
      *Knee_torque_command = t26;
    } else {
      *Knee_torque_command = -max_torque;
    }

    /* Function to prevent the desired joint angles from changing to fast.  */
    /* Works via saturation */
    if (-max_torque >= max_torque) {
      t14 = -max_torque;
    } else {
      t14 = max_torque;
    }

    if (*Ankle_torque_command <= t14) {
      t26 = *Ankle_torque_command;
    } else {
      t26 = t14;
    }

    if (t26 >= -max_torque) {
      *Ankle_torque_command = t26;
    } else {
      *Ankle_torque_command = -max_torque;
    }
  }

  /*     %% Persistent/output variable management */
  *Esys_integrate_out = 0.0F;
  knee_pos_prev = Knee_joint_position;
  ankle_pos_prev = Ankle_joint_position;
  hip_vel_prev = *hip_vel;
  hip_pos_prev = IMU_pitch;
  u_pbc_knee_prev = *U_PBC_K;
  u_pbc_ankle_prev = *U_PBC_A;
  *knee_joint_vel = knee_vel_prev;
  *ankle_joint_vel = ankle_vel_prev;
  *foot_contact = FootContact;
  *stance = Stance;
  *swing = Swing;
  *IMU_LIVE_OUT = IMU_LIVE;
  *PushOffOut = 0.0F;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void SLIP_KPBC_init(void)
{
  knee_pos_prev = 0.0F;
  ankle_pos_prev = 0.0F;
  knee_vel_prev = 0.0F;
  ankle_vel_prev = 0.0F;
  hip_vel_prev = 0.0F;
  u_pbc_knee_prev = 0.0F;
  u_pbc_ankle_prev = 0.0F;
  Stance = true;
  Swing = false;
  ForceCount = 0.0F;
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
 * File trailer for SLIP_KPBC.c
 *
 * [EOF]
 */
