#include <stdint.h>

#ifndef LINALG_FUNCTIONS_H
#define LINALG_FUNCTIONS_H

uint32_t func_solve(void *f, double *stackpos);
uint32_t func_eigvals(void *f, double *stackpos);
uint32_t func_det(void *f, double *stackpos);

#endif
