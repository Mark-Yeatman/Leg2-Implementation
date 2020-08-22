/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: PDControlTest.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 03-Oct-2019 15:42:50
 */

#ifndef PDCONTROLTEST_H
#define PDCONTROLTEST_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "PDControlTest_types.h"

/* Function Declarations */
extern void PDControlTest(float knee_pos, float ankle_pos, float knee_des, float
  ankle_des, float dt, float kp_k, float kd_k, float kp_a, float kd_a, float
  time_in, float filter_coeff, float *u_k, float *u_a, float *t_out, float
  *dt_out);
extern void knee_pos_prev_not_empty_init(void);
extern void t_prev_not_empty_init(void);

#endif

/*
 * File trailer for PDControlTest.h
 *
 * [EOF]
 */
