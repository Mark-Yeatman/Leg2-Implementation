/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xnrm2.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 25-Aug-2020 12:53:36
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC.h"
#include "xnrm2.h"

/* Function Definitions */

/*
 * Arguments    : const float x[4]
 * Return Type  : float
 */
float xnrm2(const float x[4])
{
  float y;
  float scale;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  scale = 1.17549435E-38F;
  for (k = 0; k < 2; k++) {
    absxk = (float)fabs(x[k]);
    if (absxk > scale) {
      t = scale / absxk;
      y = 1.0F + y * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * (float)sqrt(y);
}

/*
 * File trailer for xnrm2.c
 *
 * [EOF]
 */
