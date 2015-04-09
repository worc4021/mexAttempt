#pragma once
#include <gmp.h>


struct GMPmat{
    mpq_t *data;
    size_t m, n;
};

struct dMat{
    double *data;
    size_t m, n;
};