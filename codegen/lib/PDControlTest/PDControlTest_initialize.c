/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: PDControlTest_initialize.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 03-Oct-2019 15:42:50
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PDControlTest.h"
#include "PDControlTest_initialize.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void PDControlTest_initialize(void)
{
  rt_InitInfAndNaN(8U);
  t_prev_not_empty_init();
  knee_pos_prev_not_empty_init();
}

/*
 * File trailer for PDControlTest_initialize.c
 *
 * [EOF]
 */
