#pragma once
#include "GMPtypes.h"

struct dMat *readMXArray(const mxArray *pm);
mxArray *writeMXArray(const struct dMat *A);