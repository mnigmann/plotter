#include <stdint.h>

#ifndef INTERVALS_H
#define INTERVALS_H

uint32_t interval_value(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_add(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_multiply(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_div(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_exponentiate(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_equals(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_greater(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_compare_single(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_compare(void *f, double *hstackpos, double *lstackpos);

uint32_t interval_sine(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_cosine(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_tangent(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_arctan(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_exp(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_factorial(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_conjugate(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_abs(void *f, double *hstackpos, double *lstackpos);

uint32_t interval_point(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_list(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_extract_x(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_extract_y(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_user_defined(void *f, double *hstackpos, double *lstackpos);
uint32_t interval_integrate_gsl(void *f, double *hstackpos, double *lstackpos);
#endif

