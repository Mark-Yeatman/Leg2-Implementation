/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: SLIP_KPBC_v4_initialize.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 26-Sep-2019 12:16:46
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC_v4.h"
#include "SLIP_KPBC_v4_initialize.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void SLIP_KPBC_v4_initialize(void)
{
  rt_InitInfAndNaN(8U);
  t_prev_not_empty_init();
  knee_pos_prev_not_empty_init();
}

/*
 * File trailer for SLIP_KPBC_v4_initialize.c
 *
 * [EOF]
 */
