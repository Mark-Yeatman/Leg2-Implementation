/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: SLIP_rtwutil.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 15-Oct-2019 15:51:04
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP.h"
#include "SLIP_rtwutil.h"
#include "sind.h"

/* Function Definitions */

/*
 * Arguments    : float u0
 *                float u1
 * Return Type  : float
 */
float rt_remf_snf(float u0, float u1)
{
  float y;
  float b_u1;
  float tr;
  if (!((!rtIsNaNF(u0)) && (!rtIsInfF(u0)) && ((!rtIsNaNF(u1)) && (!rtIsInfF(u1)))))
  {
    y = ((real32_T)rtNaN);
  } else {
    if (u1 < 0.0F) {
      b_u1 = (float)ceil(u1);
    } else {
      b_u1 = (float)floor(u1);
    }

    if ((u1 != 0.0F) && (u1 != b_u1)) {
      tr = u0 / u1;
      if ((float)fabs(tr - rt_roundf_snf(tr)) <= FLT_EPSILON * (float)fabs(tr))
      {
        y = 0.0F;
      } else {
        y = (float)fmod(u0, u1);
      }
    } else {
      y = (float)fmod(u0, u1);
    }
  }

  return y;
}

/*
 * Arguments    : float u
 * Return Type  : float
 */
float rt_roundf_snf(float u)
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
 * File trailer for SLIP_rtwutil.c
 *
 * [EOF]
 */
