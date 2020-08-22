/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: sind.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 26-Sep-2019 12:16:46
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC_v4.h"
#include "sind.h"
#include "SLIP_KPBC_v4_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : float *x
 * Return Type  : void
 */
void sind(float *x)
{
  float absx;
  signed char n;
  if (!((!rtIsInfF(*x)) && (!rtIsNaNF(*x)))) {
    *x = ((real32_T)rtNaN);
  } else {
    *x = rt_remf_snf(*x, 360.0F);
    absx = (float)fabs(*x);
    if (absx > 180.0F) {
      if (*x > 0.0F) {
        *x -= 360.0F;
      } else {
        *x += 360.0F;
      }

      absx = (float)fabs(*x);
    }

    if (absx <= 45.0F) {
      *x *= 0.0174532924F;
      n = 0;
    } else if (absx <= 135.0F) {
      if (*x > 0.0F) {
        *x = 0.0174532924F * (*x - 90.0F);
        n = 1;
      } else {
        *x = 0.0174532924F * (*x + 90.0F);
        n = -1;
      }
    } else if (*x > 0.0F) {
      *x = 0.0174532924F * (*x - 180.0F);
      n = 2;
    } else {
      *x = 0.0174532924F * (*x + 180.0F);
      n = -2;
    }

    if (n == 0) {
      *x = (float)sin(*x);
    } else if (n == 1) {
      *x = (float)cos(*x);
    } else if (n == -1) {
      *x = -(float)cos(*x);
    } else {
      *x = -(float)sin(*x);
    }
  }
}

/*
 * File trailer for sind.c
 *
 * [EOF]
 */
