/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xscal.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 25-Aug-2020 12:53:36
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC.h"
#include "xscal.h"

/* Function Definitions */

/*
 * Arguments    : float a
 *                float x[4]
 *                int ix0
 * Return Type  : void
 */
void xscal(float a, float x[4], int ix0)
{
  int k;
  for (k = ix0; k <= ix0 + 1; k++) {
    x[k - 1] *= a;
  }
}

/*
 * File trailer for xscal.c
 *
 * [EOF]
 */
