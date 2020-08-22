/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: SLIP_initialize.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 15-Oct-2019 15:51:04
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP.h"
#include "SLIP_initialize.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void SLIP_initialize(void)
{
  rt_InitInfAndNaN(8U);
  t_prev_not_empty_init();
  knee_pos_prev_not_empty_init();
}

/*
 * File trailer for SLIP_initialize.c
 *
 * [EOF]
 */
