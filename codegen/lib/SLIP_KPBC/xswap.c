/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xswap.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 25-Aug-2020 12:53:36
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC.h"
#include "xswap.h"

/* Function Definitions */

/*
 * Arguments    : float x[4]
 * Return Type  : void
 */
void xswap(float x[4])
{
  int ix;
  int iy;
  int k;
  float temp;
  ix = 0;
  iy = 2;
  for (k = 0; k < 2; k++) {
    temp = x[ix];
    x[ix] = x[iy];
    x[iy] = temp;
    ix++;
    iy++;
  }
}

/*
 * File trailer for xswap.c
 *
 * [EOF]
 */
