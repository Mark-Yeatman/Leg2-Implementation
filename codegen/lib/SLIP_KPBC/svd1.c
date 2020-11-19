/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: svd1.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 25-Aug-2020 12:53:36
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "SLIP_KPBC.h"
#include "svd1.h"
#include "xscal.h"
#include "xrot.h"
#include "xrotg.h"
#include "xswap.h"
#include "xaxpy.h"
#include "xdotc.h"
#include "xnrm2.h"

/* Function Definitions */

/*
 * Arguments    : const float A[4]
 *                float U[4]
 *                float S[4]
 *                float V[4]
 * Return Type  : void
 */
void svd(const float A[4], float U[4], float S[4], float V[4])
{
  float b_A[4];
  float Vf[4];
  int ixstart;
  boolean_T apply_transform;
  float nrm;
  float s[2];
  int m;
  int kase;
  float ztest;
  float e[2];
  int q;
  int iter;
  float snorm;
  float rt;
  boolean_T exitg3;
  int qs;
  boolean_T exitg2;
  float f;
  float varargin_1[5];
  float mtmp;
  boolean_T exitg1;
  float sqds;
  for (ixstart = 0; ixstart < 4; ixstart++) {
    b_A[ixstart] = A[ixstart];
    Vf[ixstart] = 0.0F;
  }

  apply_transform = false;
  nrm = xnrm2(A);
  if (nrm > 0.0F) {
    apply_transform = true;
    if (A[0] < 0.0F) {
      nrm = -nrm;
    }

    if ((float)fabs(nrm) >= 9.86076132E-32F) {
      ztest = 1.0F / nrm;
      for (kase = 0; kase < 2; kase++) {
        b_A[kase] *= ztest;
      }
    } else {
      for (kase = 0; kase < 2; kase++) {
        b_A[kase] /= nrm;
      }
    }

    b_A[0]++;
    s[0] = -nrm;
  } else {
    s[0] = 0.0F;
  }

  if (apply_transform) {
    xaxpy(-(xdotc(b_A, b_A) / b_A[0]), b_A);
  }

  m = 2;
  s[1] = b_A[3];
  e[0] = b_A[2];
  e[1] = 0.0F;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    U[ixstart] = b_A[ixstart];
    U[2 + ixstart] = 0.0F;
  }

  U[3] = 1.0F;
  if (s[0] != 0.0F) {
    xaxpy(-(xdotc(U, U) / U[0]), U);
    for (ixstart = 0; ixstart < 2; ixstart++) {
      U[ixstart] = -U[ixstart];
    }

    U[0]++;
  } else {
    for (ixstart = 0; ixstart < 2; ixstart++) {
      U[ixstart] = 0.0F;
    }

    U[0] = 1.0F;
  }

  for (q = 1; q >= 0; q += -1) {
    for (ixstart = 0; ixstart < 2; ixstart++) {
      Vf[ixstart + (q << 1)] = 0.0F;
    }

    Vf[q + (q << 1)] = 1.0F;
  }

  nrm = b_A[2];
  for (q = 0; q < 2; q++) {
    if (s[q] != 0.0F) {
      rt = (float)fabs(s[q]);
      ztest = s[q] / rt;
      s[q] = rt;
      if (q + 1 < 2) {
        nrm /= ztest;
      }

      xscal(ztest, U, 1 + (q << 1));
    }

    if ((q + 1 < 2) && (nrm != 0.0F)) {
      rt = (float)fabs(nrm);
      ztest = rt / nrm;
      nrm = rt;
      s[1] *= ztest;
      xscal(ztest, Vf, 3);
    }

    e[0] = nrm;
  }

  iter = 0;
  snorm = 0.0F;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    nrm = (float)fabs(s[ixstart]);
    ztest = (float)fabs(e[ixstart]);
    if ((nrm >= ztest) || rtIsNaNF(ztest)) {
    } else {
      nrm = ztest;
    }

    if (!((snorm >= nrm) || rtIsNaNF(nrm))) {
      snorm = nrm;
    }
  }

  while ((m > 0) && (!(iter >= 75))) {
    q = m - 1;
    exitg3 = false;
    while (!(exitg3 || (q == 0))) {
      nrm = (float)fabs(e[0]);
      if ((nrm <= 1.1920929E-7F * ((float)fabs(s[0]) + (float)fabs(s[1]))) ||
          (nrm <= 9.86076132E-32F) || ((iter > 20) && (nrm <= 1.1920929E-7F *
            snorm))) {
        e[0] = 0.0F;
        exitg3 = true;
      } else {
        q = 0;
      }
    }

    if (q == m - 1) {
      kase = 4;
    } else {
      qs = m;
      ixstart = m;
      exitg2 = false;
      while ((!exitg2) && (ixstart >= q)) {
        qs = ixstart;
        if (ixstart == q) {
          exitg2 = true;
        } else {
          nrm = 0.0F;
          if (ixstart < m) {
            nrm = (float)fabs(e[0]);
          }

          if (ixstart > q + 1) {
            nrm += (float)fabs(e[0]);
          }

          ztest = (float)fabs(s[ixstart - 1]);
          if ((ztest <= 1.1920929E-7F * nrm) || (ztest <= 9.86076132E-32F)) {
            s[ixstart - 1] = 0.0F;
            exitg2 = true;
          } else {
            ixstart--;
          }
        }
      }

      if (qs == q) {
        kase = 3;
      } else if (qs == m) {
        kase = 1;
      } else {
        kase = 2;
        q = qs;
      }
    }

    switch (kase) {
     case 1:
      f = e[0];
      e[0] = 0.0F;
      kase = m - 1;
      while (kase >= q + 1) {
        xrotg(&s[0], &f, &nrm, &ztest);
        xrot(Vf, 1, 1 + ((m - 1) << 1), nrm, ztest);
        kase = 0;
      }
      break;

     case 2:
      f = e[q - 1];
      e[q - 1] = 0.0F;
      for (kase = q; kase + 1 <= m; kase++) {
        xrotg(&s[kase], &f, &nrm, &ztest);
        f = -ztest * e[kase];
        e[kase] *= nrm;
        xrot(U, 1 + (kase << 1), 1 + ((q - 1) << 1), nrm, ztest);
      }
      break;

     case 3:
      varargin_1[0] = (float)fabs(s[m - 1]);
      varargin_1[1] = (float)fabs(s[0]);
      varargin_1[2] = (float)fabs(e[0]);
      varargin_1[3] = (float)fabs(s[q]);
      varargin_1[4] = (float)fabs(e[q]);
      ixstart = 1;
      mtmp = varargin_1[0];
      if (rtIsNaNF(varargin_1[0])) {
        kase = 2;
        exitg1 = false;
        while ((!exitg1) && (kase < 6)) {
          ixstart = kase;
          if (!rtIsNaNF(varargin_1[kase - 1])) {
            mtmp = varargin_1[kase - 1];
            exitg1 = true;
          } else {
            kase++;
          }
        }
      }

      if (ixstart < 5) {
        while (ixstart + 1 < 6) {
          if (varargin_1[ixstart] > mtmp) {
            mtmp = varargin_1[ixstart];
          }

          ixstart++;
        }
      }

      f = s[m - 1] / mtmp;
      nrm = s[0] / mtmp;
      ztest = e[0] / mtmp;
      sqds = s[q] / mtmp;
      rt = ((nrm + f) * (nrm - f) + ztest * ztest) / 2.0F;
      nrm = f * ztest;
      nrm *= nrm;
      if ((rt != 0.0F) || (nrm != 0.0F)) {
        ztest = (float)sqrt(rt * rt + nrm);
        if (rt < 0.0F) {
          ztest = -ztest;
        }

        ztest = nrm / (rt + ztest);
      } else {
        ztest = 0.0F;
      }

      f = (sqds + f) * (sqds - f) + ztest;
      rt = sqds * (e[q] / mtmp);
      while (q + 1 <= 1) {
        xrotg(&f, &rt, &nrm, &ztest);
        f = nrm * s[0] + ztest * e[0];
        e[0] = nrm * e[0] - ztest * s[0];
        rt = ztest * s[1];
        s[1] *= nrm;
        xrot(Vf, 1, 3, nrm, ztest);
        s[0] = f;
        xrotg(&s[0], &rt, &nrm, &ztest);
        f = nrm * e[0] + ztest * s[1];
        s[1] = -ztest * e[0] + nrm * s[1];
        rt = ztest * e[1];
        e[1] *= nrm;
        xrot(U, 1, 3, nrm, ztest);
        q = 1;
      }

      e[0] = f;
      iter++;
      break;

     default:
      if (s[q] < 0.0F) {
        s[q] = -s[q];
        xscal(-1.0F, Vf, 1 + (q << 1));
      }

      while ((q + 1 < 2) && (s[0] < s[1])) {
        rt = s[0];
        s[0] = s[1];
        s[1] = rt;
        xswap(Vf);
        xswap(U);
        q = 1;
      }

      iter = 0;
      m--;
      break;
    }
  }

  for (kase = 0; kase < 2; kase++) {
    e[kase] = s[kase];
    for (ixstart = 0; ixstart < 2; ixstart++) {
      V[ixstart + (kase << 1)] = Vf[ixstart + (kase << 1)];
    }
  }

  for (ixstart = 0; ixstart < 4; ixstart++) {
    S[ixstart] = 0.0F;
  }

  for (kase = 0; kase < 2; kase++) {
    S[kase + (kase << 1)] = e[kase];
  }
}

/*
 * File trailer for svd1.c
 *
 * [EOF]
 */
