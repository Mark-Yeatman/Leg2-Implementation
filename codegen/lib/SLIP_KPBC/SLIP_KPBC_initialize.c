/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: SLIP_KPBC_initialize.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 01-Jul-2020 11:22:48
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC.h"
#include "SLIP_KPBC_initialize.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void SLIP_KPBC_initialize(void)
{
  rt_InitInfAndNaN(8U);
  knee_pos_prev_not_empty_init();
  SLIP_KPBC_init();
}

/*
 * File trailer for SLIP_KPBC_initialize.c
 *
 * [EOF]
 */
