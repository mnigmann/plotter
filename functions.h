#include <stdint.h>

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define TYPE_LIST 0x08
#define TYPE_DOUBLE 0x00
#define TYPE_POINT 0x01
#define TYPE_COLOR 0x02
#define TYPE_BOOLEAN 0x03
#define TYPE_ELLIPSIS 0x04
#define TYPE_MASK 0x07
#define TYPE_ABSOLUTE_ADDR 0x10

uint32_t func_div(void *f, double *stackpos);
uint32_t func_floor(void *f, double *stackpos);
uint32_t func_mod(void *f, double *stackpos);
uint32_t func_sine(void *f, double *stackpos);
uint32_t func_cosine(void *f, double *stackpos);
uint32_t func_arctan(void *f, double *stackpos);
uint32_t func_add(void *f, double *stackpos);
uint32_t func_multiply(void *f, double *stackpos);
uint32_t func_exponentiate(void *f, double *stackpos);
uint32_t func_user_defined(void *f, double *stackpos);

uint32_t func_list(void *f, double *stackpos);
uint32_t func_index(void *f, double *stackpos);
uint32_t func_point(void *f, double *stackpos);
uint32_t func_ellipsis(void *f, double *stackpos);

uint32_t func_for(void *f, double *stackpos);
uint32_t func_equals(void *f, double *stackpos);

#endif
