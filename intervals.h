#include <stdint.h>

#ifndef INTERVALS_H
#define INTERVALS_H

uint32_t interval_value(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_add(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_multiply(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_div(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_exponentiate(void *f, double *hstackpos, double *lstackpos);

uint32_t interval_sine(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_cosine(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_exp(void *f, double *hstackpos, double *lstackpos);

uint32_t interval_point(void *f, double *hstackpos, double *lstackpos);

#endif

