#include <stdint.h>
#include "parse.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define TYPE_LIST 0x08
#define TYPE_DOUBLE 0x00
#define TYPE_POINT 0x01
#define TYPE_COLOR 0x02
#define TYPE_BOOLEAN 0x03
#define TYPE_ELLIPSIS 0x04
#define TYPE_POLYGON 0x05
#define TYPE_MASK 0x07
#define TYPE_ABSOLUTE_ADDR 0x10

const oper_data *oper_lookup(uint32_t (*ptr)(void*, double*));

uint32_t func_value(void *f, double *stackpos);

uint32_t func_div(void *f, double *stackpos);
uint32_t func_floor(void *f, double *stackpos);
uint32_t func_mod(void *f, double *stackpos);
uint32_t func_max(void *f, double *stackpos);
uint32_t func_sine(void *f, double *stackpos);
uint32_t func_cosine(void *f, double *stackpos);
uint32_t func_arctan(void *f, double *stackpos);
uint32_t func_arcsin(void *f, double *stackpos);
uint32_t func_log(void *f, double *stackpos);
uint32_t func_sub(void *f, double *stackpos);
uint32_t func_exp(void *f, double *stackpos);
uint32_t func_factorial(void *f, double *stackpos);
uint32_t func_abs(void *f, double *stackpos);
uint32_t func_add(void *f, double *stackpos);
uint32_t func_multiply(void *f, double *stackpos);
uint32_t func_exponentiate(void *f, double *stackpos);
uint32_t func_user_defined(void *f, double *stackpos);

uint32_t func_list(void *f, double *stackpos);
uint32_t func_index(void *f, double *stackpos);
uint32_t func_point(void *f, double *stackpos);
uint32_t func_polygon(void *f, double *stackpos);
uint32_t func_rgb(void *f, double *stackpos);
uint32_t func_ellipsis(void *f, double *stackpos);

uint32_t func_for(void *f, double *stackpos);
uint32_t func_equals(void *f, double *stackpos);
uint32_t func_greater(void *f, double *stackpos);
uint32_t func_compare_single(void *f, double *stackpos);
uint32_t func_compare(void *f, double *stackpos);
uint32_t func_compare_sub_single(void *f, double *stackpos);
uint32_t func_compare_sub(void *f, double *stackpos);

uint32_t func_extract_x(void *f, double *stackpos);
uint32_t func_extract_y(void *f, double *stackpos);
uint32_t func_assign(void *f, double *stackpos);
uint32_t func_chain_actions(void *f, double *stackpos);
uint32_t func_sum(void *f, double *stackpos);
uint32_t func_prod(void *f, double *stackpos);
uint32_t func_total(void *f, double *stackpos);
uint32_t func_distance(void *f, double *stackpos);
uint32_t func_conditional(void *f, double *stackpos);
uint32_t func_sort(void *f, double *stackpos);
uint32_t func_join(void *f, double *stackpos);
uint32_t func_length(void *f, double *stackpos);

uint32_t func_integrate(void *f, double *stackpos);
uint32_t func_convert_polar(void *f, double *stackpos);

#endif
