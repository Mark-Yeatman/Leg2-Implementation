/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xaxpy.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 25-Aug-2020 12:53:36
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC.h"
#include "xaxpy.h"

/* Function Definitions */

/*
 * Arguments    : float a
 *                float y[4]
 * Return Type  : void
 */
void xaxpy(float a, float y[4])
{
  int ix;
  int iy;
  int k;
  if (!(a == 0.0F)) {
    ix = 0;
    iy = 2;
    for (k = 0; k < 2; k++) {
      y[iy] += a * y[ix];
      ix++;
      iy++;
    }
  }
}

/*
 * File trailer for xaxpy.c
 *
 * [EOF]
 */
