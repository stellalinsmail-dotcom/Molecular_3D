//
// File: _coder_newgetEB_api.h
//
// MATLAB Coder version            : 24.1
// C/C++ source code generated on  : 2025-10-20 19:58:52
//

#ifndef _CODER_NEWGETEB_API_H
#define _CODER_NEWGETEB_API_H

// Include Files
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <algorithm>
#include <cstring>

// Variable Declarations
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

// Function Declarations
real_T newgetEB(real_T dr_ij, real_T kb);

void newgetEB_api(const mxArray *const prhs[2], const mxArray **plhs);

void newgetEB_atexit();

void newgetEB_initialize();

void newgetEB_terminate();

void newgetEB_xil_shutdown();

void newgetEB_xil_terminate();

#endif
//
// File trailer for _coder_newgetEB_api.h
//
// [EOF]
//
