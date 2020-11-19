/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: SLIP_KPBC.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 31-Aug-2020 12:07:22
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC.h"

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
static float COP_prev;

/* Function Declarations */
static float rt_powf_snf(float u0, float u1);
static float rt_roundf_snf(float u);

/* Function Definitions */

/*
 * Arguments    : float u0
 *                float u1
 * Return Type  : float
 */
static float rt_powf_snf(float u0, float u1)
{
  float y;
  float f0;
  float f1;
  if (rtIsNaNF(u0) || rtIsNaNF(u1)) {
    y = ((real32_T)rtNaN);
  } else {
    f0 = (float)fabs(u0);
    f1 = (float)fabs(u1);
    if (rtIsInfF(u1)) {
      if (f0 == 1.0F) {
        y = ((real32_T)rtNaN);
      } else if (f0 > 1.0F) {
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
    } else if (f1 == 0.0F) {
      y = 1.0F;
    } else if (f1 == 1.0F) {
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
 * Arguments    : float u
 * Return Type  : float
 */
static float rt_roundf_snf(float u)
{
  float y;
  if ((float)fabs(u) < 8.388608E+6F) {
    if (u >= 0.5F) {
      y = (float)floor(u + 0.5F);
    } else if (u > -0.5F) {
      y = u * 0.0F;
    } else {
      y = (float)ceil(u - 0.5F);
    }
  } else {
    y = u;
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
 *                float md
 *                float Eref
 *                float knee_stop_low
 *                float knee_stop_high
 *                float ankle_stop_low
 *                float ankle_stop_high
 *                float max_torque
 *                float v0
 *                float Fric_Comp
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
 *                float *U_S_KNEE
 *                float *U_S_ANKLE
 *                float *COPFX
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
               pbc_gain_knee, float md, float Eref, float knee_stop_low, float
               knee_stop_high, float ankle_stop_low, float ankle_stop_high,
               float max_torque, float v0, float Fric_Comp, float F_thresh,
               float Command_State, float KPBC_max_torque, float Joint_Bio_Sat,
               float *Knee_torque_command, float *Ankle_torque_command, float
               *deltaL, float *hip_pos, float *Esys, float *Esys_integrate_out,
               float *U_S_KNEE, float *U_S_ANKLE, float *COPFX, float
               *U_LIN_DAMP_A, float *U_STOP_K, float *U_STOP_A, float *U_PBC_K,
               float *U_PBC_A, float *knee_des_out, float *ankle_des_out, float *
               foot_contact, float *stance, float *swing, float *phase_var_out,
               float *IMU_LIVE_OUT, float *StanceGain, float *SwingGain, float
               *knee_joint_vel, float *ankle_joint_vel, float *hip_vel, float
               *PushOffOut)
{
  float x[3];
  float M0_idx_0;
  float scale;
  int i;
  float M0_idx_2;
  boolean_T FootContact;
  float Pow_Ankle;
  float b_x[4];
  float F[3];
  float c[3];
  static const signed char a[3] = { 0, 0, 1 };

  float t6;
  float c_x;
  float t12;
  float t19;
  float J_L[2];
  boolean_T IMU_LIVE;
  signed char i0;
  float u_lin_spring[2];
  float u_lin_damp[2];
  int trueCount;
  float b_Knee_torque_command;
  float ankle_limits[2];
  static const short iv0[8] = { 1, 139, 400, 530, 720, 848, 976, 1001 };

  static const short knee_ind[8] = { 1, 139, 400, 530, 720, 848, 976, 1001 };

  float c_Knee_torque_command;
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

  float b_Ankle_torque_command;
  float c_Ankle_torque_command;
  float b_Pow_Ankle;
  float b_scale;
  float c_Pow_Ankle;
  float c_scale;

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
  /* params = [Ms,Mf,Isz,Ifz,lt,ls,lc,la,rs,rfy,rfx,COPfx,g]; */
  /*  converted lb to kg */
  /*  converted lb to kg */
  /* kg*m^2 */
  /* kg*m^2 */
  /* lt = single(0.3733); %meters */
  /* meters */
  /* meters */
  /* meters */
  /* meters */
  /* meters */
  /* meters */
  /* meters/s^2; */
  *U_LIN_DAMP_A = 0.0F;
  *U_S_KNEE = 0.0F;
  *U_S_ANKLE = 0.0F;
  *U_PBC_K = 0.0F;
  *U_PBC_A = 0.0F;

  /*  Position/velocity hard limits */
  /* deg/s */
  /* deg/s */
  /* deg/s */
  /* deg/s */
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

  /*     %% Determine foot contact */
  x[0] = Load_cell_x_force;
  x[1] = Load_cell_y_force;
  x[2] = Load_cell_z_force;
  M0_idx_0 = 0.0F;
  scale = 1.17549435E-38F;
  for (i = 0; i < 3; i++) {
    M0_idx_2 = (float)fabs(x[i]);
    if (M0_idx_2 > scale) {
      Pow_Ankle = scale / M0_idx_2;
      M0_idx_0 = 1.0F + M0_idx_0 * Pow_Ankle * Pow_Ankle;
      scale = M0_idx_2;
    } else {
      Pow_Ankle = M0_idx_2 / scale;
      M0_idx_0 += Pow_Ankle * Pow_Ankle;
    }
  }

  M0_idx_0 = scale * (float)sqrt(M0_idx_0);
  if (M0_idx_0 > F_thresh) {
    ForceCount++;
    FootContact = (ForceCount > 2.0F);
  } else {
    FootContact = false;
    ForceCount = 0.0F;
  }

  *StanceGain = 0.0F;
  *SwingGain = 0.0F;

  /* Human Leg as a Robot states */
  b_x[0] = Knee_joint_position;

  /* MAKE SURE TO MANAGE BIOMECHANICS VS RIGHT-HAND-RULE AXIS ORIENTATION */
  b_x[1] = Ankle_joint_position;
  b_x[2] = knee_vel_prev;
  b_x[3] = ankle_vel_prev;

  /* Notes: */
  /*  Knee axis treated as origin. */
  /*  Orientation of knee and ankle rotation is biomechanical */
  /* variables: */
  /*  x - 4x1 array [knee pos, ankle pos, knee vel, ankle vel] */
  /*  lt -thigh length  */
  /*  ls -shank length */
  /*  lfx, lfy - foot CoM x,y */
  /*  lc - load cell y dist from ankle axis */
  /*  la - ankle axis y dist to bottom of foot */
  /*  COPfx - COP distance along foot from ankle axis projection, to be */
  /*    computed from load cell */
  /*  Ms - shank mass */
  /*  Mf - foot mass */
  /*  COPfx calc */
  /* from pilot walking trials */
  F[0] = Load_cell_x_force;
  F[1] = Load_cell_y_force;
  F[2] = Load_cell_z_force;
  M0_idx_0 = (-0.0F * Load_cell_z_force - -0.0628F * Load_cell_y_force) +
    Load_cell_x_moment;
  scale = (-0.0628F * Load_cell_x_force - -0.0F * Load_cell_z_force) +
    Load_cell_y_moment;
  M0_idx_2 = (-0.0F * Load_cell_y_force - -0.0F * Load_cell_x_force) +
    Load_cell_z_moment;
  c[0] = 0.0F * M0_idx_2 - scale;
  c[1] = M0_idx_0 - 0.0F * M0_idx_2;
  c[2] = 0.0F * scale - 0.0F * M0_idx_0;
  Pow_Ankle = 0.0F;
  for (i = 0; i < 3; i++) {
    Pow_Ankle += F[i] * (float)a[i];
  }

  for (i = 0; i < 3; i++) {
    F[i] = c[i] / Pow_Ankle;
  }

  M0_idx_0 = c[0] / Pow_Ankle;
  if (!FootContact) {
    /* (-2,8) inches */
    M0_idx_0 = -0.0508F;
    COP_prev = -0.0508F;
  } else {
    /* only allow the COP to progress forwards */
    if (F[0] < COP_prev) {
      M0_idx_0 = COP_prev;
    } else {
      COP_prev = F[0];
    }
  }

  /* Function to prevent the desired joint angles from changing to fast.  */
  /* Works via saturation */
  if (M0_idx_0 <= 0.06F) {
    Pow_Ankle = M0_idx_0;
  } else {
    Pow_Ankle = 0.06F;
  }

  if (Pow_Ankle >= -0.0508F) {
    *COPFX = Pow_Ankle;
  } else {
    *COPFX = -0.0508F;
  }

  /* meters */
  /*      M = M_func(x,params); */
  /*      C = C_func(x,params); */
  /*      G = G_func(x,params);     */
  /* L_FUNC */
  /*     L = L_FUNC(IN1,IN2) */
  /*     This function was generated by the Symbolic Math Toolbox version 8.5. */
  /*     25-Aug-2020 13:50:39 */
  M0_idx_0 = Knee_joint_position * 3.14159274F / 180.0F;
  t6 = (Knee_joint_position + Ankle_joint_position) * 3.14159274F / 180.0F;
  c_x = (float)cos(t6);
  M0_idx_2 = (float)sin(t6);
  Pow_Ankle = ((lt - *COPFX * M0_idx_2) + 0.0628F * c_x) + (0.3733F + -lt) *
    (float)cos(M0_idx_0);
  scale = (*COPFX * c_x + 0.0628F * M0_idx_2) + (0.3733F + -lt) * (float)sin
    (M0_idx_0);
  c_x = (float)sqrt(Pow_Ankle * Pow_Ankle + scale * scale);

  /* J_L_FUNC */
  /*     JACOBL = J_L_FUNC(IN1,IN2) */
  /*     This function was generated by the Symbolic Math Toolbox version 8.5. */
  /*     25-Aug-2020 13:50:39 */
  M0_idx_0 = Knee_joint_position * 3.14159274F / 180.0F;
  t6 = Ankle_joint_position * 3.14159274F / 180.0F;
  scale = (Knee_joint_position + Ankle_joint_position) * 3.14159274F / 180.0F;
  M0_idx_2 = (float)cos(scale);
  Pow_Ankle = (float)sin(scale);
  t12 = *COPFX * M0_idx_2;
  t19 = (t12 + 0.0628F * Pow_Ankle) + (0.3733F + -lt) * (float)sin(M0_idx_0);
  scale = ((lt + 0.0628F * M0_idx_2) + -(*COPFX * Pow_Ankle)) + (0.3733F + -lt) *
    (float)cos(M0_idx_0);
  scale = 1.0F / (float)sqrt(t19 * t19 + scale * scale);
  J_L[0] = -lt * t19 * scale;
  J_L[1] = -scale * ((0.0628F * (lt * Pow_Ankle + (0.3733F + -lt) * (float)sin
    (t6)) + lt * t12) + *COPFX * (0.3733F + -lt) * (float)cos(t6));

  /* J_dot_L = J_dot_L_func(x,params);     */
  /*      Theta = Theta_func(x,params); */
  /*      J_Theta = J_Theta_func(x,params); */
  /*      J_dot_Theta = J_dot_Theta_func(x,params);    */
  /*      J_z = [J_L;J_Theta]; */
  /*      J_z_dot = [J_dot_L;J_dot_Theta];       */
  /*      J_C = J_C_func(x,params);    */
  *deltaL = c_x - L0;
  M0_idx_0 = 0.0F;
  for (i = 0; i < 2; i++) {
    M0_idx_0 += J_L[i] * b_x[2 + i];
  }

  if (md < 10.0F) {
    md = 10.0F;
  } else {
    if (md > 10000.0F) {
      md = 10000.0F;
    }
  }

  /* M_z = diag([md,md]); */
  /* Calculate system energy */
  Pow_Ankle = c_x - L0;

  /* virtual spring potential energy */
  /* virtual mass kinetic energy */
  *Esys = 0.5F * k * (Pow_Ankle * Pow_Ankle) + 0.5F * md * (M0_idx_0 * M0_idx_0);

  /* Phase Variable  */
  /* s_a=clamp(1+(1-s_m)/(q_h_0-q_h_m)*(hip-q_h_0),0,1); */
  Pow_Ankle = 0.4F * (IMU_pitch - -20.0F) / 43.0F;
  if (0.6F + Pow_Ankle <= 1.0F) {
    Pow_Ankle += 0.6F;
  } else {
    Pow_Ankle = 1.0F;
  }

  if (Pow_Ankle >= 0.6F) {
    *phase_var_out = Pow_Ankle;
  } else {
    *phase_var_out = 0.6F;
  }

  /* IMU State check */
  IMU_LIVE = ((float)fabs(hip_pos_prev - IMU_pitch) < 1.0F);

  /*     %% Phase State management, set autodetection or command stance/swing */
  scale = rt_roundf_snf(Command_State);
  if (scale < 128.0F) {
    if (scale >= -128.0F) {
      i0 = (signed char)scale;
    } else {
      i0 = MIN_int8_T;
    }
  } else if (scale >= 128.0F) {
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
      scale = (time_in - t_switched_phase) / 0.24000001F;
      if (scale <= 1.0F) {
        *StanceGain = scale;
      } else {
        *StanceGain = 1.0F;
      }

      *SwingGain = 1.0F - *StanceGain;
    } else {
      if (Swing) {
        scale = (time_in - t_switched_phase) / 0.24000001F;
        if (scale <= 1.0F) {
          *SwingGain = scale;
        } else {
          *SwingGain = 1.0F;
        }

        *StanceGain = 1.0F - *SwingGain;
      }
    }
  } else {
    scale = rt_roundf_snf(Command_State);
    if (scale < 128.0F) {
      if (scale >= -128.0F) {
        i0 = (signed char)scale;
      } else {
        i0 = MIN_int8_T;
      }
    } else if (scale >= 128.0F) {
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
      scale = rt_roundf_snf(Command_State);
      if (scale < 128.0F) {
        if (scale >= -128.0F) {
          i0 = (signed char)scale;
        } else {
          i0 = MIN_int8_T;
        }
      } else if (scale >= 128.0F) {
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

  /*     %% Torque computation */
  if (SLIP_ON != 0.0F) {
    /*         %% Stance                 */
    /* Virtual spring */
    /* MAKE SURE TO MANAGE BIOMECHANICS VS RIGHT-HAND-RULE AXIS ORIENTATION */
    Pow_Ankle = -k * *deltaL;
    scale = -d * M0_idx_0;
    for (i = 0; i < 2; i++) {
      u_lin_spring[i] = Pow_Ankle * J_L[i];
      u_lin_damp[i] = scale * J_L[i];
    }

    *U_LIN_DAMP_A = u_lin_damp[1];
    t12 = u_lin_damp[0] + u_lin_spring[0];
    t19 = u_lin_damp[1] + u_lin_spring[1];
    if (KPBC_ON != 0.0F) {
      /* Also use energy tracking controller    */
      Pow_Ankle = pbc_gain_knee * (*Esys - Eref) * M0_idx_0;
      for (i = 0; i < 2; i++) {
        J_L[i] *= Pow_Ankle;
      }

      /* MAKE SURE TO MANAGE BIOMECHANICS VS RIGHT-HAND-RULE AXIS ORIENTATION */
      *U_PBC_K = J_L[0];
      *U_PBC_A = J_L[1];

      /* Moving average filter and saturation laws on value and rate of torque        */
      /* Absolute saturation */
      if (KPBC_max_torque > 0.0F) {
        /* Function to prevent the desired joint angles from changing to fast.  */
        /* Works via saturation */
        if (-KPBC_max_torque >= KPBC_max_torque) {
          M0_idx_2 = -KPBC_max_torque;
        } else {
          M0_idx_2 = KPBC_max_torque;
        }

        if (J_L[0] <= M0_idx_2) {
          Pow_Ankle = J_L[0];
        } else {
          Pow_Ankle = M0_idx_2;
        }

        if (Pow_Ankle >= -KPBC_max_torque) {
          *U_PBC_K = Pow_Ankle;
        } else {
          *U_PBC_K = -KPBC_max_torque;
        }

        /* Function to prevent the desired joint angles from changing to fast.  */
        /* Works via saturation */
        if (-KPBC_max_torque >= KPBC_max_torque) {
          M0_idx_2 = -KPBC_max_torque;
        } else {
          M0_idx_2 = KPBC_max_torque;
        }

        if (J_L[1] <= M0_idx_2) {
          Pow_Ankle = J_L[1];
        } else {
          Pow_Ankle = M0_idx_2;
        }

        if (Pow_Ankle >= -KPBC_max_torque) {
          *U_PBC_A = Pow_Ankle;
        } else {
          *U_PBC_A = -KPBC_max_torque;
        }
      }

      /* Rate saturation */
      if (KPBC_max_torque_rate > 0.0F) {
        scale = (*U_PBC_K - u_pbc_knee_prev) / dt;
        M0_idx_0 = (*U_PBC_A - u_pbc_ankle_prev) / dt;

        /* Function to prevent the desired joint angles from changing to fast.  */
        /* Works via saturation */
        if (-KPBC_max_torque_rate >= KPBC_max_torque_rate) {
          M0_idx_2 = -KPBC_max_torque_rate;
        } else {
          M0_idx_2 = KPBC_max_torque_rate;
        }

        if (scale <= M0_idx_2) {
          Pow_Ankle = scale;
        } else {
          Pow_Ankle = M0_idx_2;
        }

        /* Function to prevent the desired joint angles from changing to fast.  */
        /* Works via saturation */
        if (-KPBC_max_torque_rate >= KPBC_max_torque_rate) {
          M0_idx_2 = -KPBC_max_torque_rate;
        } else {
          M0_idx_2 = KPBC_max_torque_rate;
        }

        if (M0_idx_0 <= M0_idx_2) {
          scale = M0_idx_0;
        } else {
          scale = M0_idx_2;
        }

        if (Pow_Ankle >= -KPBC_max_torque_rate) {
          c_Pow_Ankle = Pow_Ankle;
        } else {
          c_Pow_Ankle = -KPBC_max_torque_rate;
        }

        *U_PBC_K = c_Pow_Ankle * dt + u_pbc_knee_prev;
        if (scale >= -KPBC_max_torque_rate) {
          c_scale = scale;
        } else {
          c_scale = -KPBC_max_torque_rate;
        }

        *U_PBC_A = c_scale * dt + u_pbc_ankle_prev;
      }

      /* Exponential smoothing of torque */
      *U_PBC_K = (1.0F - KPBC_filter_coeff) * u_pbc_knee_prev +
        KPBC_filter_coeff * *U_PBC_K;
      *U_PBC_A = (1.0F - KPBC_filter_coeff) * u_pbc_ankle_prev +
        KPBC_filter_coeff * *U_PBC_A;

      /* Biomemetic power saturation */
      if (Joint_Bio_Sat != 0.0F) {
        scale = knee_vel_prev * *U_PBC_K;
        Pow_Ankle = ankle_vel_prev * *U_PBC_A;
        if (Pow_Ankle < 0.0F) {
          b_Pow_Ankle = -1.0F;
        } else if (Pow_Ankle > 0.0F) {
          b_Pow_Ankle = 1.0F;
        } else if (Pow_Ankle == 0.0F) {
          b_Pow_Ankle = 0.0F;
        } else {
          b_Pow_Ankle = Pow_Ankle;
        }

        if (b_Pow_Ankle == -1.0F) {
          *U_PBC_A = 0.0F;
        }

        if (scale < 0.0F) {
          b_scale = -1.0F;
        } else if (scale > 0.0F) {
          b_scale = 1.0F;
        } else if (scale == 0.0F) {
          b_scale = 0.0F;
        } else {
          b_scale = scale;
        }

        if (b_scale == 1.0F) {
          *U_PBC_K = 0.0F;
        }
      }

      t12 = (t12 + *U_PBC_K) - u_lin_damp[0];
      t19 = (t19 + *U_PBC_A) - u_lin_damp[1];
    }

    /*         %% Swing */
    if (IMU_LIVE) {
      /* Use hip pos as phase variable to index winters data */
      Pow_Ankle = *phase_var_out * 1000.0F;
      trueCount = 0;
      for (i = 0; i < 8; i++) {
        if (Pow_Ankle - (float)knee_ind[i] >= 0.0F) {
          trueCount++;
        }
      }

      c_x = (Pow_Ankle - (float)iv0[(int)(float)trueCount - 1]) / (float)(iv0
        [(int)((float)trueCount + 1.0F) - 1] - iv0[(int)(float)trueCount - 1]);
      *knee_des_out = ((a_knee[(int)(float)trueCount - 1] * ((c_x - 1.0F) * (c_x
        - 1.0F)) + a_knee[(int)(float)trueCount + 6] * rt_powf_snf(c_x - 1.0F,
        6.0F)) + a_knee[(int)(float)trueCount + 20] * (c_x - 1.0F)) + a_knee
        [(int)(float)trueCount + 13];
      trueCount = 0;
      for (i = 0; i < 11; i++) {
        if (Pow_Ankle - (float)ankle_ind[i] >= 0.0F) {
          trueCount++;
        }
      }

      c_x = (Pow_Ankle - (float)iv1[(int)(float)trueCount - 1]) / (float)(iv1
        [(int)((float)trueCount + 1.0F) - 1] - iv1[(int)(float)trueCount - 1]);
      *ankle_des_out = ((a_ankle[(int)(float)trueCount - 1] * ((c_x - 1.0F) *
        (c_x - 1.0F)) + a_ankle[(int)(float)trueCount + 9] * rt_powf_snf(c_x -
        1.0F, 6.0F)) + a_ankle[(int)(float)trueCount + 29] * (c_x - 1.0F)) +
        a_ankle[(int)(float)trueCount + 19];

      /* follow the trjactory via PD controller */
      /* This control basically assumes that the desired joint angles are static i.e., the desired velocities are all zero. */
      /*  Helper functions */
      scale = kp_knee * (*knee_des_out - Knee_joint_position) + kd_knee *
        -knee_vel_prev;
      M0_idx_0 = kp_ankle * (*ankle_des_out - Ankle_joint_position) + kd_ankle *
        -ankle_vel_prev;
    } else {
      *knee_des_out = Knee_joint_position;
      *ankle_des_out = Ankle_joint_position;

      /* This control basically assumes that the desired joint angles are static i.e., the desired velocities are all zero. */
      /*  Helper functions */
      scale = kp_knee * (Knee_joint_position - Knee_joint_position) + kd_knee *
        -knee_vel_prev;
      M0_idx_0 = kp_ankle * (Ankle_joint_position - Ankle_joint_position) +
        kd_ankle * -ankle_vel_prev;
    }

    /* Smooth between stance and swing torques */
    *Knee_torque_command = *StanceGain * t12 + *SwingGain * scale;
    *Ankle_torque_command = *StanceGain * t19 + *SwingGain * M0_idx_0;
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

  /*     %% Friction Compensation */
  /* v0 = 0.01; */
  if (Fric_Comp != 0.0F) {
    if (knee_vel_prev < 0.0F) {
      scale = 0.0F;
    } else if (knee_vel_prev < v0) {
      scale = knee_vel_prev / v0;
    } else {
      scale = 1.0F;
    }

    if (knee_vel_prev < -v0) {
      M0_idx_0 = 0.0F;
    } else if (knee_vel_prev < 0.0F) {
      M0_idx_0 = (knee_vel_prev - (-v0)) / (0.0F - (-v0));
    } else {
      M0_idx_0 = 1.0F;
    }

    if (*Knee_torque_command < 0.0F) {
      b_Knee_torque_command = 0.0F;
    } else if (*Knee_torque_command < 5.0F) {
      b_Knee_torque_command = *Knee_torque_command / 5.0F;
    } else {
      b_Knee_torque_command = 1.0F;
    }

    if (*Knee_torque_command < -5.0F) {
      c_Knee_torque_command = 0.0F;
    } else if (*Knee_torque_command < 0.0F) {
      c_Knee_torque_command = (*Knee_torque_command - -5.0F) / 5.0F;
    } else {
      c_Knee_torque_command = 1.0F;
    }

    *Knee_torque_command += (scale * 0.8F + (1.0F - M0_idx_0) * -1.7F) + (1.0F -
      scale) * M0_idx_0 * (b_Knee_torque_command * 1.6F + (1.0F -
      c_Knee_torque_command) * -2.8F);
    if (ankle_vel_prev < 0.0F) {
      scale = 0.0F;
    } else if (ankle_vel_prev < v0) {
      scale = ankle_vel_prev / v0;
    } else {
      scale = 1.0F;
    }

    if (ankle_vel_prev < -v0) {
      M0_idx_0 = 0.0F;
    } else if (ankle_vel_prev < 0.0F) {
      M0_idx_0 = (ankle_vel_prev - (-v0)) / (0.0F - (-v0));
    } else {
      M0_idx_0 = 1.0F;
    }

    if (*Ankle_torque_command < 0.0F) {
      b_Ankle_torque_command = 0.0F;
    } else if (*Ankle_torque_command < 5.0F) {
      b_Ankle_torque_command = *Ankle_torque_command / 5.0F;
    } else {
      b_Ankle_torque_command = 1.0F;
    }

    if (*Ankle_torque_command < -5.0F) {
      c_Ankle_torque_command = 0.0F;
    } else if (*Ankle_torque_command < 0.0F) {
      c_Ankle_torque_command = (*Ankle_torque_command - -5.0F) / 5.0F;
    } else {
      c_Ankle_torque_command = 1.0F;
    }

    *Ankle_torque_command += (scale * 0.8F + (1.0F - M0_idx_0) * -1.7F) + (1.0F
      - scale) * M0_idx_0 * (b_Ankle_torque_command * 1.6F + (1.0F -
      c_Ankle_torque_command) * -2.8F);
  }

  /*     %% Virtual hard stops */
  for (i = 0; i < 2; i++) {
    u_lin_spring[i] = 0.0F;
  }

  if ((knee_stop_low < knee_stop_high) && (ankle_stop_low < ankle_stop_high)) {
    J_L[0] = knee_stop_low;
    J_L[1] = knee_stop_high;

    /* deg, need to make sure its ordered */
    ankle_limits[0] = ankle_stop_low;
    ankle_limits[1] = ankle_stop_high;
  } else {
    /* deg, need to make sure its ordered */
    for (i = 0; i < 2; i++) {
      J_L[i] = 2.0F + 103.0F * (float)i;
      ankle_limits[i] = -35.0F + 70.0F * (float)i;
    }
  }

  if (Knee_joint_position < J_L[0]) {
    u_lin_spring[0] = -kp_knee * (Knee_joint_position - J_L[0]) - kd_knee *
      knee_vel_prev;
  }

  if (J_L[1] < Knee_joint_position) {
    u_lin_spring[0] = -kp_knee * (Knee_joint_position - J_L[1]) - kd_knee *
      knee_vel_prev;
  }

  if (Ankle_joint_position < ankle_limits[0]) {
    u_lin_spring[1] = -kp_ankle * (Ankle_joint_position - ankle_limits[0]) -
      kd_ankle * ankle_vel_prev;
  }

  if (ankle_limits[1] < Ankle_joint_position) {
    u_lin_spring[1] = -kp_ankle * (Ankle_joint_position - ankle_limits[1]) -
      kd_ankle * ankle_vel_prev;
  }

  *U_STOP_K = u_lin_spring[0];
  *U_STOP_A = u_lin_spring[1];
  *Knee_torque_command += u_lin_spring[0];
  *Ankle_torque_command += u_lin_spring[1];

  /*     %% Saturate torque output */
  if (max_torque > 0.0F) {
    /* Function to prevent the desired joint angles from changing to fast.  */
    /* Works via saturation */
    if (-max_torque >= max_torque) {
      M0_idx_2 = -max_torque;
    } else {
      M0_idx_2 = max_torque;
    }

    if (*Knee_torque_command <= M0_idx_2) {
      Pow_Ankle = *Knee_torque_command;
    } else {
      Pow_Ankle = M0_idx_2;
    }

    if (Pow_Ankle >= -max_torque) {
      *Knee_torque_command = Pow_Ankle;
    } else {
      *Knee_torque_command = -max_torque;
    }

    /* Function to prevent the desired joint angles from changing to fast.  */
    /* Works via saturation */
    if (-max_torque >= max_torque) {
      M0_idx_2 = -max_torque;
    } else {
      M0_idx_2 = max_torque;
    }

    if (*Ankle_torque_command <= M0_idx_2) {
      Pow_Ankle = *Ankle_torque_command;
    } else {
      Pow_Ankle = M0_idx_2;
    }

    if (Pow_Ankle >= -max_torque) {
      *Ankle_torque_command = Pow_Ankle;
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
  COP_prev = -0.0508F;
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
