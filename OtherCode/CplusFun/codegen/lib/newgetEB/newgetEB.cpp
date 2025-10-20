//
// File: newgetEB.cpp
//
// MATLAB Coder version            : 24.1
// C/C++ source code generated on  : 2025-10-20 19:58:52
//

// Include Files
#include "newgetEB.h"

// Function Definitions
//
// Arguments    : double dr_ij
//                double kb
// Return Type  : double
//
double newgetEB(double dr_ij, double kb)
{
  double EB_tmp;
  EB_tmp = dr_ij * dr_ij;
  return 143.9525 * kb / 2.0 * EB_tmp *
         ((-2.0 * dr_ij + 1.0) + 2.3333333333333335 * EB_tmp);
}

//
// File trailer for newgetEB.cpp
//
// [EOF]
//
