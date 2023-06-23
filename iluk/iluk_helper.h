#ifndef ILUK_HELPER
#define ILUK_HELPER

#include "mex.h"

void sparse_sqr_transpose(double*, mwIndex*, mwIndex*, mwSize);
void reallocNzmax(mxArray*, double**, mwIndex**, mwSize);

#endif