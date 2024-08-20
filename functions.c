#include "functions.h"
#include "parse.h"
#include "config.h"
#include "intervals.h"
#include "linalg_functions.h"
#include "linalg_intervals.h"
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#define SIGN_BIT(v) (((v->value_type)&0x80) ? -1 : 1)
#define VALUE(a) ((((a->value_type)&0x40) ? *(((variable*)(a->value))->pointer) : *((double*)(a->value))) * SIGN_BIT(a))
#define VALUE_LIST(a, idx) ((((a->value_type)&0x40) ? ((double*)(((variable*)(a->value))->pointer))[idx] : *((double*)(a->value) + idx)) * SIGN_BIT(a))

#define FAIL(...) {printf(__VA_ARGS__); exit(EXIT_FAILURE);}
#define GET_STEP(type) (step_table[(type) & TYPE_MASK])
#define IS_TYPE(type, ref) (((type) & TYPE_MASK) == (ref))

const static double lanczos_table[9] = {
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

const static uint32_t step_table[8] = {1, 2, 3, 1, 1, MAX_POLYGON_SIZE*2, 0, 0};


uint32_t func_general_one_arg(void *f, double *stackpos, double (*oper)(double)) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t len, type;
    type = arg->oper(arg, stackpos);
    len = type>>8;
    for (int i=0; i < len; i++) stackpos[i] = oper(stackpos[i])*SIGN_BIT(fs);
    return (len<<8) | (type & TYPE_LIST);
}

uint32_t func_general_two_args(void *f, double *stackpos, double (*oper)(double, double)) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t t1, t2, l1, l2;
    double *val1, *val2;
    t1 = arg->oper(arg, stackpos);
    val1 = stackpos;
    l1 = t1>>8;
    arg = arg->next_arg;
    val2 = stackpos+(l1?l1:1);
    t2 = arg->oper(arg, val2);
    l2 = t2>>8;
    //printf("func_general_two_args: "); print_object(t1, val1); printf(" and "); print_object(t2, val2); printf("\n");
    if (!(t1 & TYPE_LIST) && !(t2 & TYPE_LIST)) {
        //printf("func_general_two_args: %p (%f), %p (%f)", val1, *val1, val2, *val2);
        stackpos[0] = oper(*val1, *val2)*SIGN_BIT(fs);
        //printf(" --> %f\n", stackpos[0]);
        return 1<<8;
    } else if ((t1 & TYPE_LIST) && !(t2 & TYPE_LIST)) {
        double v2 = *val2;
        for (int i=0; i < l1; i++) stackpos[i] = oper(val1[i], v2)*SIGN_BIT(fs);
        return t1;
    } else if (!(t1 & TYPE_LIST) && (t2 & TYPE_LIST)) {
        double v1 = *val1;
        for (int i=0; i < l2; i++) stackpos[i] = oper(v1, val2[i])*SIGN_BIT(fs);
        return t2;
    } else if (l1 > l2) {
        for (int i=0; i < l2; i++) stackpos[i] = oper(val1[i], val2[i])*SIGN_BIT(fs);
        return t2;
    } else {
        for (int i=0; i < l1; i++) stackpos[i] = oper(val1[i], val2[i])*SIGN_BIT(fs);
        return t1;
    }
}

uint32_t func_value(void *f, double *stackpos) {
    function *fs = (function*)f; 
    uint32_t type;
    double *ptr;
    variable *var;
    if (fs->value_type&0x40) {
        var = ((variable*)(fs->value));
        type = var->type;
        ptr = var->pointer;
        if (!ptr) FAIL("ERROR: function block %p references variable %p (%s) with null pointer\n", fs, fs->value, var->name);
        // If the variable is an action, run its source code
        if (var->flags & VARIABLE_ACTION) {
            printf("Action found: %p, %p\n", var, var->pointer);
            ((function*)(var->pointer))->oper(var->pointer, stackpos);
        }
    } else {
        type = fs->value_type;
        ptr = (double*)(fs->value);
        if (!ptr) FAIL("ERROR: pointer in function block %p is null\n", fs);
    }
    uint32_t len = type>>8;
    for (int i=0; i < len; i++) stackpos[i] = ptr[i]*SIGN_BIT(fs);
    return type;
}

double mdiv(double n, double d) {
    return n/d;
}

double mmod(double n, double d) {
    double i = fmod(n, d);
    if (n && ((n>0) ^ (d>0))) i += d;
    return i;
}

double mmax(double a, double b) {
    if (a>b) return a;
    return b;
}

double mmin(double a, double b) {
    if (a<b) return a;
    return b;
}

double mlogfac(double x) {
    double v = lanczos_table[0];
    for (int i=1; i < 9; i++) {
        v += lanczos_table[i]/(x+i);
    }
    // log(2*pi)/2
    return 0.9189385332046727 + (x+0.5)*log(x+7.5) - (x+7.5) + log(v);
}

double mfac(double x) {
    if (x < -0.5) return -M_PI / sin(M_PI*x) * exp(-mlogfac(-1-x));
    return exp(mlogfac(x));
}

double msub(double x, double y) {
    return x-y;
}

double mround2(double x, double n) {
    return round(x*pow(10, n))/pow(10, n);
}

double msign(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

uint32_t func_sub(void *f, double *stackpos) {
    return func_general_two_args(f, stackpos, msub);
}

uint32_t func_div(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t t1, t2, l1, l2;
    double *val1, *val2;
    uint8_t step1, step2;
    t1 = arg->oper(arg, stackpos);
    val1 = stackpos;
    l1 = t1>>8;
    step1 = GET_STEP(t1);

    arg = arg->next_arg;
    val2 = stackpos+l1;
    t2 = arg->oper(arg, val2);
    l2 = t2>>8;
    step2 = GET_STEP(t2);

    uint32_t result_len = l1/step1;
    double val;
    if (((l2/step2 < result_len) || !(t1 & TYPE_LIST)) && (t2 & TYPE_LIST)) result_len = l2/step2;
    if (IS_TYPE(t1, TYPE_DOUBLE)) {
        if (IS_TYPE(t2, TYPE_DOUBLE)) {
            // Dividing real numbers
            if (!(t1 & TYPE_LIST)) {
                val = val1[0];
                val1 = &val;
            }
            for (int i=0; i < result_len; i++) stackpos[i] = val1[i%l1]/val2[i%l2]*SIGN_BIT(fs);
            return (result_len<<8) | ((t1 | t2) & 0xff);
        } else {
            // divide real by complex
            for (int i=l1+l2-1; i>=0; i--) stackpos[2*result_len+i]=stackpos[i];
            val1 += 2*result_len;
            val2 += 2*result_len;
            complex double zn;
            for (int i=0; i < result_len; i++) {
                zn = val1[i%l1];
                zn /= (val2[(2*i)%l2] + val2[(2*i+1)%l2]*I);
                stackpos[2*i] = creal(zn)*SIGN_BIT(fs);
                stackpos[2*i+1] = cimag(zn)*SIGN_BIT(fs);
            }
            return (result_len<<9) | ((t1 | t2) & 0xff);
        }
    } else {
        for (int i=l1+l2-1; i>=0; i--) stackpos[2*result_len+i]=stackpos[i];
        val1 += 2*result_len;
        val2 += 2*result_len;
        if (IS_TYPE(t2, TYPE_DOUBLE)) {
            // divide complex by real
            //printf("dividing complex by real: "); print_object(t1, val1); printf(" / "); print_object(t2, val2); printf("\n");
            complex double zn;
            for (int i=0; i < result_len; i++) {
                zn = val1[(2*i)%l1] + I*val1[(2*i+1)%l1];
                zn /= val2[i%l2];
                stackpos[2*i] = creal(zn)*SIGN_BIT(fs);
                stackpos[2*i+1] = cimag(zn)*SIGN_BIT(fs);
            }
        } else {
            complex double zn;
            for (int i=0; i < result_len; i++) {
                zn = val1[(2*i)%l1] + I*val1[(2*i+1)%l1];
                zn /= val2[(2*i)%l2] + I*val2[(2*i+1)%l2];
                stackpos[2*i] = creal(zn)*SIGN_BIT(fs);
                stackpos[2*i+1] = cimag(zn)*SIGN_BIT(fs);
            }
        }
        return (result_len<<9) | ((t1 | t2) & 0xff);
    }
}

uint32_t func_mod(void *f, double *stackpos) {
    return func_general_two_args(f, stackpos, mmod);
}

uint32_t func_floor(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, floor);
}

uint32_t func_round(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    if (arg->next_arg) return func_general_two_args(f, stackpos, mround2);
    else return func_general_one_arg(f, stackpos, round);
}

uint32_t func_sine(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, sin);
}

uint32_t func_cosine(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, cos);
}

uint32_t func_tangent(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, tan);
}

uint32_t func_factorial(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, mfac);
}

uint32_t func_abs(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t type = arg->oper(arg, stackpos);
    uint32_t len = type>>8;
    if (IS_TYPE(type, TYPE_POINT)) {
        for (int i=0, j=0; i < len; i+=2, j++) stackpos[j] = hypot(stackpos[i], stackpos[i+1]);
        type = ((len/2)<<8) | (type & TYPE_LIST);
        return type;
    } else {
        for (int i=0; i<len; i++) {
            if (stackpos[i] < 0) stackpos[i] = -stackpos[i];
        }
        return type;
    }
}

uint32_t func_log(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    // logarithms of factorials must be computed separately for better performance
    if (arg->oper == func_factorial) {
        uint32_t type, len;
        type = arg->first_arg->oper(arg->first_arg, stackpos);
        len = type>>8;
        for (int i=0; i < len; i++) stackpos[i] = SIGN_BIT(fs)*mlogfac(stackpos[i]);
        return type;
    }
    return func_general_one_arg(f, stackpos, log);
}

uint32_t func_exp(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, exp);
}

uint32_t func_max(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    if (arg->next_arg) {
        uint32_t type = func_general_two_args(f, stackpos, mmax);
        return type;
    }
    
    uint32_t type = arg->oper(arg, stackpos);
    for (int i=1; i < (type>>8); i++) {
        if (stackpos[i] > stackpos[0]) stackpos[0] = stackpos[i];
    }
    stackpos[0] *= SIGN_BIT(fs);
    return (1<<8);
}

uint32_t func_min(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    if (arg->next_arg) {
        uint32_t type = func_general_two_args(f, stackpos, mmin);
        return type;
    }
    
    uint32_t type = arg->oper(arg, stackpos);
    for (int i=1; i < (type>>8); i++) {
        if (stackpos[i] < stackpos[0]) stackpos[0] = stackpos[i];
    }
    stackpos[0] *= SIGN_BIT(fs);
    return (1<<8);
}

uint32_t func_sign(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, msign);
}

// The arctan function can have either one or two arguments
uint32_t func_arctan(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    if (arg->next_arg) return func_general_two_args(f, stackpos, atan2);
    else return func_general_one_arg(f, stackpos, atan);
}

uint32_t func_arcsin(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, asin);
}

uint32_t func_conjugate(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;

    uint32_t type = arg->oper(arg, stackpos);
    uint32_t len = type>>8;
    if (IS_TYPE(type, TYPE_POINT)) {
        for (uint32_t i=1; i < len; i+=2) stackpos[i] *= -1;
    }
    for (uint32_t i=0; i < len; i++) stackpos[i] *= SIGN_BIT(fs);
    return type;
}

uint8_t add_in_place(double *acc, double *val, uint32_t *result_type, uint32_t *result_length, uint32_t type) {
    uint32_t len = type>>8;
    if ((*result_type & TYPE_MASK) != (type & TYPE_MASK)) {
        // Error
        return 1;
    }
    //printf("Adding (%08x): ", type); print_object(type, val); printf(" to "); print_object(*result_type, acc); printf("\n");
    if (!(*result_type & TYPE_LIST) && (type & TYPE_LIST)) {
        // To add a list to a scalar, we swap the operands and then copy down.
        uint8_t step = GET_STEP(*result_type);
        for (int i=0; i < len; i++) val[i] += acc[i%step];
        *result_length = len;
        *result_type = type;
        for (int i=0; i < len; i++) acc[i] = val[i];
        //printf("result "); print_object(*result_type, acc); printf("\n");
    } else {
        // Otherwise, we broadcast and add
        if ((type & TYPE_LIST)  && (len < *result_length)) *result_length = len;
        for (int i=0; i < *result_length; i++) acc[i] += val[i%len];
    }
    return 0;
}

uint32_t func_add(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    double sum = 0;
    uint32_t type;
    uint32_t len;
    uint32_t result_type = 0;
    uint32_t result_length = 0;
    int i=0;
    result_type = arg->oper(arg, stackpos);
    result_length = result_type>>8;
    arg = arg->next_arg;
    double *pos;
    //printf("func_add\n");
    while (arg) {
        type = arg->oper(arg, stackpos+result_length);
        //printf("adding %08x: ", type); print_object(type, stackpos+result_length); printf("\n");
        // Keep track of the position in case result_length is changed.
        pos = stackpos+result_length;
        if (add_in_place(stackpos, pos, &result_type, &result_length, type)) {
            printf("ERROR: only objects with the same type may be added (function block %p)\n", fs);
            printf("    result: "); print_object((result_length<<8) | (result_type & 0xff), stackpos);
            printf("\n    expr:   "); print_object(type, pos);
            FAIL("\n");
        }
        arg = arg->next_arg;
    }
    for (i=0; i < result_length; i++) stackpos[i] *= SIGN_BIT(fs);
    
    return (result_length<<8) | (result_type&0xff);
}

void multiply_in_place(double *stackpos, double *pos, uint32_t *result_type, uint32_t *result_length, uint32_t type) {
    uint32_t len = type>>8;
    if ((!(*result_type & TYPE_LIST) || (len/GET_STEP(type) < *result_length)) && (type & TYPE_LIST)) *result_length = len/GET_STEP(type);
    uint8_t step = GET_STEP(*result_type);
    if (IS_TYPE(*result_type, TYPE_POINT) && IS_TYPE(type, TYPE_POINT)) {
        if (!IS_TYPE(*result_type, TYPE_LIST) && IS_TYPE(type, TYPE_LIST)) {
            // Result is scalar, must convert
            complex double zn;
            complex double zr = stackpos[0] + stackpos[1]*I;
            for (int i=0; i < len; i+=2) {
                zn = stackpos[2+i] + stackpos[3+i]*I;
                zn *= zr;
                stackpos[i] = creal(zn);
                stackpos[i+1] = cimag(zn);
            }
        } else {
            complex double zn;
            len = len / 2;
            for (int i=0; i<*result_length; i++) {
                zn = stackpos[2*i] + stackpos[2*i+1]*I;
                zn *= (pos[2*(i%len)] + pos[2*(i%len) + 1]*I);
                stackpos[2*i] = creal(zn);
                stackpos[2*i+1] = cimag(zn);
            }
        }
    }
    //printf("Multiplying (%08x): ", type); print_object(type, stackpos+result_length); printf("\n");
    else if (!(*result_type & TYPE_LIST) && (type & TYPE_LIST)) {
        // To multiply a scalar/point by a list, we multiply the list by
        // the scalar/point in-place and then copy the list down. If the 
        // existing value is a point, then the list must be a list of scalars
        // and the final result will be the Kronecker product of the list
        // and the point: (x, y) * [a, b, c] = [(ax, ay), (bx, by), (cx, cy)]
        for (int i=len-1; i >= 0; i--) {
            for (int j=step-1; j >= 0; j--) pos[i*step+j] = pos[i]*stackpos[j];
        }
        *result_length = len;
        *result_type = type | *result_type;
        for (int i=0; i < len*step; i++) stackpos[i] = pos[i];
    } else {
        // Otherwise, we broadcast and add. We also know that if type is a
        // list, then result_type is as well.
        if ((type & TYPE_LIST) && (len/step < *result_length)) *result_length = len/step;
        if (IS_TYPE(type, TYPE_POINT)) {
            // If type is a point, then result_type is not a point
            if (type & TYPE_LIST) {
                // List of scalars times a list of points. result_length has already
                // been shortened to the appropriate length.
                for (int i=0; i < *result_length; i++) {
                    pos[2*i] *= stackpos[i];
                    pos[2*i+1] *= stackpos[i];
                }
                memcpy(stackpos, pos, (*result_length)*step*sizeof(double));
            } else {
                // Scalar or list of scalars times a single point
                double px = pos[0], py = pos[1];
                for (int i=*result_length-1; i>=0; i--) {
                    // y comes first because x might overwrite stackpos[0]
                    stackpos[2*i+1] = stackpos[i]*py;
                    stackpos[2*i] = stackpos[i]*px;
                }
            }
            *result_type |= TYPE_POINT;
        } else {
            for (int i=0; i < *result_length; i++) {
                for (int j=0; j < step; j++) stackpos[i*step+j] *= pos[i%len];
            }
        }
    }
}

uint32_t func_multiply(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t type;
    uint32_t len;
    uint32_t result_type = 0;
    uint32_t result_length = 0;
    int i=0;
    result_type = arg->oper(arg, stackpos);
    result_length = result_type>>8;
    result_length /= GET_STEP(result_type);
    arg = arg->next_arg;
    //printf("Initial value (%08x): ", result_type); print_object(result_type, stackpos); printf("\n");
    double *pos;
    while (arg) {
        type = arg->oper(arg, stackpos+result_length*GET_STEP(result_type));
        //printf("multiplying %08x ", type); print_object(type, stackpos+result_length*GET_STEP(result_type)); printf("\n");
        // Keep track of the position in case result_length is changed.
        pos = stackpos+result_length*GET_STEP(result_type);
        multiply_in_place(stackpos, pos, &result_type, &result_length, type);
        arg = arg->next_arg;
    }
    for (i=0; i < result_length; i++) stackpos[i] *= SIGN_BIT(fs);
    
    return ((GET_STEP(result_type)*result_length)<<8) | (result_type&0xff);
}

uint32_t func_exponentiate(void *f, double *stackpos) {
    return func_general_two_args(f, stackpos, pow);
}

uint32_t func_user_defined(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    function *target = fs->value;
    // If target is an operation, then target->next_arg will have been assigned to point to the variable for the function's first argument
    // As many variables will be filled as were given to the function
    variable *target_arg = (variable*)(target->next_arg);
    uint32_t st = 0;
    uint32_t type, len;
    uint32_t varnum = 0;
    //printf("user_defined, target %p, target arg %p\n", target, target_arg);
    do {
        //printf("argument points to function block %p, value_type %08x, value %p\n", arg, arg->value_type, arg->value);
        if (arg->oper) {
            type = arg->oper(arg, stackpos+st);
            target_arg[varnum].pointer = stackpos+st;
            target_arg[varnum].type = type;
            st += (type>>8);
        } else {
            type = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
            len = type>>8;
            for (int i=0; i < len; i++) stackpos[st+i] = VALUE_LIST(arg, i);
            target_arg[varnum].pointer = stackpos+st;
            target_arg[varnum].type = type;
            st += len;
        }
        //printf("argument is %08x: ", type); print_object(type, target_arg[varnum].pointer); printf("\n");
        varnum++;
        arg = arg->next_arg;
    } while (arg);
    type = target->oper(target, stackpos+st);
    len = type>>8;
    // Shrink the stack by removing unused values
    for (int i=0; i < len; i++) stackpos[i] = stackpos[st+i]*SIGN_BIT(fs);
    //printf("return type is %08x: ", type); print_object(type, stackpos); printf("\n");
    return type;
}

uint32_t func_list(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    if (!arg) return TYPE_LIST;
    uint8_t type = 0xff;
    uint32_t argtype;
    uint32_t arglen;
    uint32_t st=0;
    uint8_t ellipsis = 0;
    do {
        argtype = arg->oper(arg, stackpos+st);
        if ((argtype == TYPE_ELLIPSIS) && ((type & TYPE_MASK) == TYPE_DOUBLE) && (arg->next_arg)) {
            ellipsis = 2;
        } else {
            if (argtype & TYPE_LIST) FAIL("ERROR: Cannot store a list inside a list\n");
            if ((type != 0xff) && ((argtype&0xff) != type)) FAIL("ERROR: All list elements must have the same type\n");
            st += argtype>>8;
            type = argtype & 0xff;
        }
        if (ellipsis == 2) ellipsis--;
        else if (ellipsis) {
            ellipsis = 0;
            if (stackpos[st-1] != stackpos[st-2]) {
                double end = stackpos[st-1];
                double step = (st > 2 ? stackpos[st-2]-stackpos[st-3] : (stackpos[st-2] > end ? -1 : 1));
                if (step && ((step > 0) ^ (stackpos[st-1] < stackpos[st-2]))) {
                    st--;
                    for (double val=stackpos[st-1]; val < end; st++) {
                        val += step;
                        stackpos[st] = val;
                    }
                }
            }
            // If the start and end values are equal, delete one so ony one of them
            // occurs in the result. Thus, \left[1,...,1\right] evaluates to [1]
            else st--;
        }
        arg = arg->next_arg;
    } while (arg);
    for (int i=0; i < st; i++) stackpos[i] *= SIGN_BIT(fs);
    return (st<<8) | type | TYPE_LIST;
}


uint32_t func_index(void *f, double *stackpos) {
    // Three ways of indexing a list:
    //  * Boolean indexing, where the list is truncated as necesary -> list
    //  * By list of indices -> list
    //  * Singe numerical index -> value
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t t1, t2, l1, l2;
    double *val1, *val2;
    t1 = arg->oper(arg, stackpos);
    l1 = t1 >> 8;
    val1 = stackpos;
    if (!(t1 & TYPE_LIST)) {
        print_object(t1, val1); printf("\n");
        FAIL("ERROR: Cannot index a non-list (first argument %p)\n", arg);
    }
    
    arg = arg->next_arg;
    t2 = arg->oper(arg, stackpos+l1);
    l2 = t2>>8;
    val2 = stackpos + l1;

    uint8_t step = GET_STEP(t1);
    //printf("indexing "); print_object(t1, val1); printf(" with (%08x) ", t2); print_object(t2, val2); printf("\n");
    if ((t2 & TYPE_LIST) && ((t2 & TYPE_MASK) == TYPE_DOUBLE)) {
        // List of indices
        int32_t v2;
        for (int i=0; i < l2; i++) {
            v2 = step*((int32_t)val2[i]-1);
            for (uint8_t j=0; j < step; j++) stackpos[l1+l2+i*step+j] = (((v2 >= l1) || (v2 < 0)) ? NAN : val1[v2+j]);
        }
        for (int i=0; i < l2*step; i++) stackpos[i] = stackpos[l1+l2+i]*SIGN_BIT(fs);
        return ((step*l2)<<8) | (t1 & 0xff);
    } else if ((t2 & TYPE_LIST) && ((t2 & TYPE_MASK) == TYPE_BOOLEAN)) {
        // List of booleans
        uint32_t lmin = (l1 < l2 ? l1 : l2);
        uint32_t lr = 0;
        for (int i=0; i < lmin; i++) {
            if (val2[i] == 0) continue;
            for (uint8_t j=0; j < step; j++) stackpos[lr+j] = val1[i*step+j]*SIGN_BIT(fs);
            lr+=step;
        }
        return (lr << 8) | (t1 & 0xff);
    } else if ((t2 & TYPE_MASK) == TYPE_DOUBLE) {
        int32_t v2 = step*(val2[0]-1);
        for (uint8_t j=0; j < step; j++) stackpos[j] = (((v2 >= l1) || (v2 < 0)) ? NAN : val1[v2+j])*SIGN_BIT(fs);
        return (step << 8) | (t1 & 0xff & ~TYPE_LIST);
    } else FAIL("ERROR: invalid indexing operation: type %08x indexes type %08x\n", t2, t1);
}

uint32_t func_point(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t t1, t2, l1, l2;
    double *val1, *val2;
    if (arg->oper) {
        t1 = arg->oper(arg, stackpos);
        l1 = t1 >> 8;
    } else {
        t1 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
        l1 = t1 >> 8;
        for (int i=0; i < l1; i++) stackpos[i] = VALUE_LIST(arg, i);
    }
    val1 = stackpos;
    arg = arg->next_arg;
    if (arg->oper) {
        t2 = arg->oper(arg, stackpos+l1);
        l2 = t2>>8;
    } else {
        t2 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
        l2 = t2>>8;
        for (int i=0; i < l2; i++) stackpos[l1+i] = VALUE_LIST(arg, i);
    }
    val2 = stackpos + l1;
    if (!(t1 & TYPE_LIST) && !(t2 & TYPE_LIST)) {
        stackpos[0] *= SIGN_BIT(fs);
        stackpos[1] *= SIGN_BIT(fs);
        return (2<<8) | TYPE_POINT;
    } else if ((t1 & TYPE_LIST) && !(t2 & TYPE_LIST)) {
        double v2 = *val2;
        for (int i=l1-1; i >= 0; i--) {
            stackpos[2*i] = stackpos[i]*SIGN_BIT(fs);
            stackpos[2*i+1] = v2;
        }
        return (l1<<9) | TYPE_POINT | TYPE_LIST;
    } else if (!(t1 & TYPE_LIST) && (t2 & TYPE_LIST)) {
        double v1 = *val1;
        for (int i=l2-1; i >= 0; i--) {
            stackpos[2*i] = v1;
            stackpos[2*i+1] = stackpos[i+1]*SIGN_BIT(fs);
        }
        return (l2<<9) | TYPE_POINT | TYPE_LIST;
    } else if (l1 > l2) {
        for (int i=0; i < l2; i++) stackpos[l1+l2+i] = stackpos[i];
        val1 = stackpos+l1+l2;
        for (int i=0; i < l2; i++) {
            stackpos[2*i] = val1[i]*SIGN_BIT(fs);
            stackpos[2*i+1] = val2[i]*SIGN_BIT(fs);
        }
        return (l2<<9) | TYPE_POINT | TYPE_LIST;
    } else {
        // Copy l1, which may be shorter, onto the end of the stack
        for (int i=0; i < l1; i++) stackpos[l1+l2+i] = stackpos[i];
        val1 = stackpos+l1+l2;
        for (int i=0; i < l1; i++) {
            stackpos[2*i] = val1[i]*SIGN_BIT(fs);
            stackpos[2*i+1] = val2[i]*SIGN_BIT(fs);
        }
        return (l1<<9) | TYPE_POINT | TYPE_LIST;
    }
    
}

uint32_t func_polygon(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    
    int n_args = 0;
    uint32_t arglen, argtype;
    while (arg) {
        n_args++;
        arg = arg->next_arg;
    }
    arg = fs->first_arg;

    if (n_args == 1) {
        argtype = arg->oper(arg, stackpos);
        arglen = argtype >> 8;
        if (arglen > MAX_POLYGON_SIZE) FAIL("ERROR: polygon has too many points\n");
        return (arglen<<8) | TYPE_POLYGON;
    }
    
    if (n_args > MAX_POLYGON_SIZE) FAIL("ERROR: polygon has too many points\n");
    uint32_t first_len=0, st=0, min_len=0;
    while (arg) {
        argtype = arg->oper(arg, stackpos+st);
        if ((argtype & TYPE_MASK) != TYPE_POINT) FAIL("ERROR: polygon arguments must be points");
        arglen = argtype>>8;
        //printf("Points are: "); print_object(argtype, stackpos+st); printf("\n");
        if (first_len == 0) {
            first_len = arglen;
            min_len = arglen;
        } else if (arglen < min_len) min_len = arglen;
        st += first_len;
        arg = arg->next_arg;
    }
    // Lazy, slow solution
    double *temp = malloc(first_len*MAX_POLYGON_SIZE*sizeof(double));
    memcpy(temp, stackpos, first_len*MAX_POLYGON_SIZE*sizeof(double));
    for (int i=0; i < (first_len>>1); i++) {
        for (int j=0; j < n_args; j++) {
            stackpos[i*2*MAX_POLYGON_SIZE + 2*j] = temp[j*first_len + 2*i];
            stackpos[i*2*MAX_POLYGON_SIZE + 2*j + 1] = temp[j*first_len + 2*i + 1];
        }
        for (int j=n_args; j < MAX_POLYGON_SIZE; j++) {
            stackpos[i*2*MAX_POLYGON_SIZE + 2*j] = NAN;
            stackpos[i*2*MAX_POLYGON_SIZE + 2*j + 1] = NAN;
        }
    }
    free(temp);
    return ((min_len*MAX_POLYGON_SIZE)<<8) | TYPE_LIST | TYPE_POLYGON;
}

uint32_t func_rgb(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    
    uint32_t r, g, b;
    uint32_t rl, gl, bl;
    uint32_t result_len = 0;
    uint32_t result_type = 0;
    uint32_t st=0, argtype, arglen;
    argtype = arg->oper(arg, stackpos);
    rl = argtype>>8;
    result_len = rl;
    st += rl;

    arg = arg->next_arg;
    g = st;
    argtype = arg->oper(arg, stackpos+st);
    gl = argtype>>8;
    if ((argtype & TYPE_LIST) && ((result_type & TYPE_LIST) ? (gl < result_len) : 1)) {
        result_len = gl;
        result_type = argtype;
    }
    st += gl;

    arg = arg->next_arg;
    b = st;
    argtype = arg->oper(arg, stackpos+st);
    bl = argtype>>8;
    if ((argtype & TYPE_LIST) && ((result_type & TYPE_LIST) ? (bl < result_len) : 1)) {
        result_len = bl;
        result_type = argtype;
    }
    st += bl;

    double *temp = malloc(st*sizeof(double));
    memcpy(temp, stackpos, st*sizeof(double));
    for (int i=0; i < result_len; i++) {
        stackpos[3*i] = temp[i%rl];
        stackpos[3*i+1] = temp[g + i%gl];
        stackpos[3*i+2] = temp[b + i%bl];
    }
    free(temp);
    return ((3*result_len)<<8) | (result_type & TYPE_LIST) | TYPE_COLOR;
}

uint32_t func_ellipsis(void *f, double *stackpos) {
    return TYPE_ELLIPSIS;
}

uint32_t func_for(void *f, double *stackpos) {
    /* For loops may have several variables. To evaluate a for loop, the following steps must take place
     *  * Determine the number of variables
     *  * Allocate space at the beginning of the stack for ending addresses for each list
     *  * Evaluate all of the lists to be iterated over and store them on the stack
     *  * Iterate until all of the pointers are at their ends
     *  * Copy the result into position.
     * 
     */
    function *fs = (function*)f;
    function *arg1 = fs->first_arg;
    function *arg = arg1->next_arg;
    uint8_t type = 0xff;
    uint32_t argtype;
    uint32_t arglen;
    uint8_t ellipsis = 0;
    uint8_t n_args = 0;
    while (arg) {
        n_args++;
        arg = arg->next_arg;
    }
    uint32_t st=2*n_args;
    uint32_t *types;
    variable *first_var = NULL;

    arg = arg1->next_arg;
    function *var_value;
    int i = 0;
    uint32_t total = 1;
    do {
        if (first_var == NULL) first_var = arg->first_arg->value;
        var_value = arg->first_arg->next_arg;
        if (var_value->oper) {
            argtype = var_value->oper(var_value, stackpos+st);
        } else {
            argtype = ((var_value->value_type)&0x40 ? ((variable*)(var_value->value))->type : var_value->value_type);
            arglen = argtype>>8;
            for (int j=0; j < arglen; j++) stackpos[st+j] = VALUE_LIST(var_value, j);
        }
        first_var[i].pointer = stackpos+st;
        types = (uint32_t*)(stackpos+2*i);
        *types = argtype>>8;
        types = (uint32_t*)(stackpos+2*i+1);
        *types = 1;
        if (argtype & TYPE_POINT) *types = 2;
        first_var[i].type = ((*types)<<8) | (argtype & 0xff & ~TYPE_LIST);
        st += (argtype>>8);
        if (argtype & TYPE_POINT) total *= (argtype >> 9);
        else total *= (argtype >> 8);
        arg = arg->next_arg;
        i++;
    } while (arg);

    if ((arg1->oper == NULL) && !(arg1->value_type & 0x40)) {
        st = 0;
        double *val = arg1->value;
        if (arg1->value_type & TYPE_POINT) {
            for (i=0; i < total; i++) {
                stackpos[2*i] = val[0];
                stackpos[2*i+1] = val[1];
            }
            return (total<<9) | (arg1->value_type & TYPE_MASK) | TYPE_LIST;
        }
        for (i=0; i < total; i++) {
            stackpos[i] = val[0];
        }
        return (total<<8) | TYPE_LIST;
    }
    
    int input_len = st;
    for (i=0; i < total; i++) {
        //printf("evaluating for expression at %f, function pointer %p\n", first_var[0].pointer[0], arg1);
        argtype = arg1->oper(arg1, stackpos+st);
        st += (argtype>>8);
        double *pos = stackpos + n_args*2;
        for (int j=0; j < n_args; j++) {
            // Step the pointer
            types = (uint32_t*)(stackpos+2*j+1);
            first_var[j].pointer+=*types;
            // Determine the end of the variable's list
            types = (uint32_t*)(stackpos+2*j);
            pos += *types;
            if (first_var[j].pointer >= pos) first_var[j].pointer -= *types;
            else break;
        }
    }
    for (i=0; i < n_args; i++) first_var[i].pointer = 0;
    for (i=input_len; i < st; i++) stackpos[i-input_len] = stackpos[i]*SIGN_BIT(fs);
    
    return ((st - input_len)<<8) | (0xff&argtype) | TYPE_LIST;
}

double equals(double a, double b) {
    return (a==b ? 1.0 : 0.0);
}

double mgt(double a, double b) {
    return (a>b ? 1.0 : 0.0);
}

double mge(double a, double b) {
    return (a>=b ? 1.0 : 0.0);
}

double mlt(double a, double b) {
    return (a<b ? 1.0 : 0.0);
}

double mle(double a, double b) {
    return (a<=b ? 1.0 : 0.0);
}

uint32_t func_equals(void *f, double *stackpos) {
    //printf("checking equals\n");
    uint32_t res = func_general_two_args(f, stackpos, equals);
    return res | TYPE_BOOLEAN;
}

uint32_t func_greater(void *f, double *stackpos) {
    //printf("checking equals\n");
    uint32_t res = func_general_two_args(f, stackpos, mgt);
    return res | TYPE_BOOLEAN;
}

uint32_t func_compare_single(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint64_t mp = *((uint64_t*)(fs->value));
    if (mp == 0) return func_general_two_args(f, stackpos, mgt) | TYPE_BOOLEAN;
    if (mp == 1) return func_general_two_args(f, stackpos, mlt) | TYPE_BOOLEAN;
    if (mp == 2) return func_general_two_args(f, stackpos, mge) | TYPE_BOOLEAN;
    if (mp == 3) return func_general_two_args(f, stackpos, mle) | TYPE_BOOLEAN;
    return 0;
}

uint32_t func_compare_sub_single(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint64_t mp = *((uint64_t*)(fs->value));
    uint32_t old_type = fs->value_type;
    // Set the sign bit if the operation is lt or le
    fs->value_type = mp<<7;
    uint32_t type = func_general_two_args(f, stackpos, msub);
    fs->value_type = old_type;
    return type;
}

uint32_t func_compare(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint64_t cmp = *((uint64_t*)(fs->value));
    double (*cmp_oper)(double, double);
    uint32_t result_type = 0;
    uint32_t result_len = 0;
    uint32_t last_len = 0;
    uint32_t type, len, newlen;
    double result=1, last;
    result_type = arg->oper(arg, stackpos);
    len = result_type>>8;
    memcpy(stackpos+len, stackpos, len*sizeof(double));
    last_len = len; result_len = len;
    for (int i=0; i < result_len; i++) stackpos[i] = 1;
    arg = arg->next_arg;
    while (arg) {
        switch (cmp & 3) {
            case 0:
                cmp_oper = mgt;
                break;
            case 1:
                cmp_oper = mlt;
                break;
            case 2:
                cmp_oper = mge;
                break;
            case 3:
                cmp_oper = mle;
                break;
        }
        cmp = cmp >> 2;
        type = arg->oper(arg, stackpos+result_len+last_len);
        len = type>>8;
        if ((type & TYPE_LIST) && !(result_type & TYPE_LIST)) {
            // If the result is not a list but the operand is, expand the result array
            for (int i=last_len+len-1; i >= 0; i--) stackpos[len+i] = stackpos[result_len+i];
            for (int i=1; i < len; i++) stackpos[i] = stackpos[0];
            result_len = len;
            result_type |= TYPE_LIST;
        }
        newlen = result_len;
        if (!(result_type & TYPE_LIST) || ((type & TYPE_LIST) && (len < result_len))) newlen = len;
        for (int i=newlen-1; i >= 0; i--) {
            stackpos[i] = stackpos[i] * cmp_oper(stackpos[result_len+i%last_len], stackpos[result_len+last_len+i%len]);
        }
        for (int i=0; i < len; i++) stackpos[newlen+i] = stackpos[result_len+last_len+i];
        result_len = newlen;
        last_len = len;
        arg = arg->next_arg;
    }
    if (result_len == 0) {
        stackpos[0] = result;
        result_len = 1;
    }
    return (result_len << 8) | (result_type & TYPE_LIST) | TYPE_BOOLEAN;
}

uint32_t func_compare_sub(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint64_t cmp = *((uint64_t*)(fs->value));
    double (*cmp_oper)(double, double);
    uint32_t result_type = 0;
    uint32_t result_len = 0;
    uint32_t last_len = 0;
    uint32_t type, len, newlen;
    double result=1, last;
    result_type = arg->oper(arg, stackpos);
    len = result_type>>8;
    memcpy(stackpos+len, stackpos, len*sizeof(double));
    last_len = len; result_len = len;
    arg = arg->next_arg;
    while (arg) {
        type = arg->oper(arg, stackpos+result_len+last_len);
        len = type>>8;
        if ((type & TYPE_LIST) && !(result_type & TYPE_LIST)) {
            // If the result is not a list but the operand is, expand the result array
            for (int i=last_len+len-1; i >= 0; i--) stackpos[len+i] = stackpos[result_len+i];
            for (int i=1; i < len; i++) stackpos[i] = stackpos[0];
            result_len = len;
            result_type |= TYPE_LIST;
        }
        newlen = result_len;
        if (!(result_type & TYPE_LIST) || ((type & TYPE_LIST) && (len < result_len))) newlen = len;
        int i2;
        if (cmp & 1) {
            // less than
            for (int i=newlen-1; i >= 0; i--) {
                i2 = result_len+i%last_len;
                stackpos[i2] = stackpos[result_len+last_len+i%len] - stackpos[i2];
                if (stackpos[i2] < stackpos[i]) stackpos[i] = stackpos[i2];
            }
        } else {
            // greater than
            for (int i=newlen-1; i >= 0; i--) {
                i2 = result_len+i%last_len;
                stackpos[i2] = stackpos[i2] - stackpos[result_len+last_len+i%len];
                if (stackpos[i2] < stackpos[i]) stackpos[i] = stackpos[i2];
            }
        }
        cmp = cmp >> 2;
        for (int i=0; i < len; i++) stackpos[newlen+i] = stackpos[result_len+last_len+i];
        result_len = newlen;
        last_len = len;
        arg = arg->next_arg;
    }
    if (result_len == 0) {
        stackpos[0] = result;
        result_len = 1;
    }
    return (result_len << 8) | (result_type & TYPE_LIST);
}

uint32_t func_extract_x(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    //printf("extracting x, arg %p, oper %p (%p), value %p, value_type %08x\n", arg, arg->oper, (func_value), arg->value, arg->value_type);
    
    uint32_t argtype, arglen;
    argtype = arg->oper(arg, stackpos);
    arglen = argtype>>8;
    //printf("extracting from type %08x\n", argtype);
    if (!((argtype & TYPE_MASK) == TYPE_POINT)) FAIL("ERROR: cannot access x-coordinate of type %08x, block %p\n", argtype, fs);
    for (int i=0; i < arglen; i += 2) {
        stackpos[i>>1] = SIGN_BIT(fs)*stackpos[i];
    }
    return ((argtype>>1) & 0xffffff00) | ((argtype&0xff) & ~TYPE_MASK);
}

uint32_t func_extract_y(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    
    uint32_t argtype, arglen;
    argtype = arg->oper(arg, stackpos);
    arglen = argtype>>8;
    if (!((argtype & TYPE_MASK) == TYPE_POINT)) FAIL("ERROR: cannot access y-coordinate of type %08x\n", argtype);
    for (int i=0; i < arglen; i += 2) {
        stackpos[i>>1] = SIGN_BIT(fs)*stackpos[i+1];
    }
    return ((argtype>>1) & 0xffffff00) | ((argtype&0xff) & ~TYPE_MASK);
}

uint32_t func_assign(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *target = fs->first_arg;
    function *value = target->next_arg;

    variable *target_var;
    if ((target->oper == func_value) && (target->value_type & 0x40)) target_var = target->value;
    else FAIL("ERROR: Cannot assign to non-variable\n");
    
    uint32_t argtype = value->oper(value, stackpos);
    double *temp = malloc((argtype>>8)*sizeof(double));
    memcpy(temp, stackpos, (argtype>>8)*sizeof(double));
    target_var->new_pointer = temp;
    target_var->new_type = argtype;
    printf("Assigning "); print_object(argtype, temp); printf(" to %s\n", target_var->name);

    return 0;
}

uint32_t func_chain_actions(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *target = fs->first_arg;
    while (target) {
        target->oper(target, stackpos);
        target = target->next_arg;
    }
    return 0;
}

uint32_t func_sum(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *var = fs->first_arg;
    function *start = var->next_arg;
    function *end = start->next_arg;
    function *expr = end->next_arg;
    variable *varptr = var->value;
    uint32_t t_start, t_end, t_expr, st=0;
    t_start = start->oper(start, stackpos);
    st += t_start>>8;
    t_end = end->oper(end, stackpos+st);
    double *val1 = stackpos, *val2 = stackpos+st;
    st += t_end>>8;
    
    uint32_t result_length = 1;
    uint32_t result_type = 0;
    stackpos[st] = 0;
    double value;
    int i;
    double start_val = *val1;
    double end_val = *val2;
    double *pos;
    uint32_t type, len;
    //printf("summing from %e to %e\n", start_val, end_val);
    if (!(t_start & TYPE_LIST) && !(t_end & TYPE_LIST)) {
        if (start_val > end_val) {
            stackpos[0] = 0;
            return 1<<8;
        }
        varptr->pointer = &value;
        varptr->type = 0x100;
        for (value=start_val; value <= end_val; value++) {
            if (value == start_val) {
                result_type = expr->oper(expr, stackpos+st);
                result_length = result_type>>8;
                continue;
            }
            type = expr->oper(expr, stackpos+result_length+st);
            // Keep track of the position in case result_length is changed.
            pos = stackpos+result_length+st;
            if (add_in_place(stackpos+st, pos, &result_type, &result_length, type)) {
                printf("ERROR: only objects with the same type may be added (function block %p)\n", fs);
                printf("    result: "); print_object((result_length<<8) | (result_type & 0xff), stackpos+st);
                printf("\n    expr:   "); print_object(type, pos);
                FAIL("\n");
            }
        }
        for (i=0; i < result_length; i++) stackpos[i] = stackpos[st+i]*SIGN_BIT(fs);
        //printf("sum is %e\n", stackpos[0]);
    }
    varptr->pointer = 0;
    return (result_length<<8) | result_type;
}

uint32_t func_prod(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *var = fs->first_arg;
    function *start = var->next_arg;
    function *end = start->next_arg;
    function *expr = end->next_arg;
    variable *varptr = var->value;
    uint32_t t_start, t_end, t_expr, st=0;
    t_start = start->oper(start, stackpos);
    st += t_start>>8;
    t_end = end->oper(end, stackpos+st);
    double *val1 = stackpos, *val2 = stackpos+st;
    st += t_end>>8;
    
    uint32_t result_length = 0;
    uint32_t result_type = 0;
    double value;
    double prod = 1;
    int i;
    double start_val = *val1;
    double end_val = *val2;
    //printf("multiplying from %e to %e\n", start_val, end_val);
    if (!(t_start & TYPE_LIST) && !(t_end & TYPE_LIST)) {
        varptr->pointer = &value;
        varptr->type = 0x100;
        if (start_val > end_val) {
            stackpos[0] = SIGN_BIT(fs);
            return 1<<8;
        }
        if (fs->value_type == (1<<8)) {
            // Special case of a scalar product of real numbers
            value = start_val;
            expr->oper(expr, stackpos+st);
            prod = stackpos[st];
            value++;
            for (; value <= end_val; value++) {
                expr->oper(expr, stackpos+st);
                prod *= stackpos[st];
            }
            stackpos[0] = prod*SIGN_BIT(fs);
            return 1<<8;
        }
        for (value=start_val; value <= end_val; value++) {
            if (value == start_val) {
                result_type = expr->oper(expr, stackpos+st);
                result_length = result_type>>8;
                continue;
            }
            t_expr = expr->oper(expr, stackpos+st+result_length);
            multiply_in_place(stackpos+st, stackpos+st+result_length, &result_type, &result_length, t_expr);
        }
    }
    varptr->pointer = 0;
    for (int j=0; j < result_length; j++) stackpos[j] = stackpos[st+j]*SIGN_BIT(fs);
    //printf("product is %e\n", stackpos[0]);
    fs->value_type = (result_length<<8) | result_type;
    printf("storing value_type %08x\n", fs->value_type);
    return (result_length<<8) | result_type;
}

uint32_t func_total(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    
    uint32_t argtype = arg->oper(arg, stackpos);
    uint32_t arglen = argtype>>8;
    if ((argtype & TYPE_MASK) == TYPE_POINT) {
        for (int i=2; i < arglen; i+=2) {
            stackpos[0] += stackpos[i];
            stackpos[1] += stackpos[i+1];
        }
        stackpos[0] *= SIGN_BIT(fs);
        stackpos[1] *= SIGN_BIT(fs);
        return (2<<8) | TYPE_POINT;
    } else {
        for (int i=1; i < arglen; i++) stackpos[0] += stackpos[i];
        stackpos[0] *= SIGN_BIT(fs);
        return (1<<8);
    }
}

uint32_t func_distance(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t t1, t2, l1, l2;
    double *val1, *val2;
    if (arg->oper) {
        t1 = arg->oper(arg, stackpos);
        val1 = stackpos;
    } else {
        t1 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
        val1 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->pointer : (double*)(arg->value));
    }
    l1 = t1>>9;
    arg = arg->next_arg;
    if (arg->oper) {
        val2 = stackpos+(l1<<1);
        t2 = arg->oper(arg, val2);
    } else {
        t2 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
        val2 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->pointer : (double*)(arg->value));
    }
    l2 = t2>>9;
    if (((t1 & TYPE_MASK) != TYPE_POINT) || ((t2 & TYPE_MASK) != TYPE_POINT)) FAIL("ERROR: distance can only be computed between points\n");
    if (!(t1 & TYPE_LIST) && !(t2 & TYPE_LIST)) {
        stackpos[0] = hypot(val1[0]-val2[0], val1[1]-val2[1])*SIGN_BIT(fs);
        return 1<<8;
    } else if ((t1 & TYPE_LIST) && !(t2 & TYPE_LIST)) {
        double v2x = val2[0], v2y = val2[1];
        for (int i=0; i < l1; i++) stackpos[i] = hypot(val1[2*i]-v2x, val1[2*i+1]-v2y)*SIGN_BIT(fs);
        return (l1<<8) | TYPE_LIST;
    } else if (!(t1 & TYPE_LIST) && (t2 & TYPE_LIST)) {
        double v1x = val1[0], v1y = val1[1];
        for (int i=0; i < l2; i++) stackpos[i] = hypot(v1x-val2[2*i], v1y-val2[2*i+1])*SIGN_BIT(fs);
        return (l2<<8) | TYPE_LIST;
    } else if (l1 > l2) {
        for (int i=0; i < l2; i++) stackpos[i] = hypot(val1[2*i]-val2[2*i], val1[2*i+1]-val2[2*i+1])*SIGN_BIT(fs);
        return (l2<<8) | TYPE_LIST;
    } else {
        for (int i=0; i < l1; i++) stackpos[i] = hypot(val1[2*i]-val2[2*i], val1[2*i+1]-val2[2*i+1])*SIGN_BIT(fs);
        return (l1<<8) | TYPE_LIST;
    }
}

uint32_t func_conditional(void *f, double *stackpos) {
    // In order to speed up computation, conditionals are short-circuited.
    // If the first conditional evaluates to true, then the conditional
    // returns 1 even if a later condition or value is a list. If all
    // values in a list are determined, evaluation stops. In normal use,
    // this should not matter, as conditional expressions generally have
    // the same type for each branch. Broadcasting occurs as needed.
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    
    double last_mask = 1;
    double result = NAN;
    double *value_ptr = NULL;
    double *result_ptr = NULL;
    double *last_mask_ptr = NULL;
    uint32_t argtype=-1, lasttype=0, lastlen=0, arglen, st=0, result_type=0, result_len=-1;
    uint8_t step, result_step = 1, skipped=0;
    //printf("Conditional found at %p, %p\n", fs, arg);
    // After first condition:
    // v- st
    // | mask_1 |
    // After first value:
    // v- st
    // | mask | value_1 | mask_1
    // | mask | result_1 |
    //                   ^- st
    // After second condition
    // | mask | result_1 | mask_2 |
    //                            ^- st+lastlen
    // After second value
    // | mask | result_1 | mask_2 | value_2 |
    // | mask | result_2
    // After third value (default)
    // | mask | result_2 | value_3 |
    // | ~mask | result_2 | value_3 |
    // | ~mask | result_3 |
    // | result_3 |

    if (!arg) {
        stackpos[0] = SIGN_BIT(fs);
        return 1<<8;
    }
    if (!(arg->next_arg)) {
        argtype = arg->oper(arg, stackpos);
        for (int i=0; i < (argtype>>8); i++) {
            if (stackpos[i]) stackpos[i] = SIGN_BIT(fs);
            else stackpos[i] = NAN;
        }
        return argtype;
    }

    while (arg) {
        skipped = 0;
        if ((argtype == ((1<<8)|TYPE_BOOLEAN)) && (last_mask_ptr[0] == 0) && (arg->oper == func_assign)) {
            // If the previous expression was a single false boolean, then we don't actually need to evaluate
            // This is useful especially when the branch of the conditional is an action
            skipped = 1;
            argtype = 1<<8;
        }
        else argtype = arg->oper(arg, stackpos+st+lastlen);
        arglen = argtype>>8;
        step = GET_STEP(argtype);
        //printf("evaluating component of conditional: ");
        //if (skipped) printf("<skipped>\n");
        //else {printf("%08x, %08x ", argtype, result_type); print_object(argtype, stackpos+st+lastlen); printf("\n");}
        // In this case, result_len refers to the number of elements in the list,
        // not the number of doubles the list takes up

        // If the result length and type have not been set, set them
        if (result_len == -1) {
            result_len = arglen;
            if (argtype & TYPE_LIST) result_type |= TYPE_LIST;
        }
        // If the existing data is not a list, then we may need to expand it
        else if ((argtype & TYPE_LIST) && !(result_type & TYPE_LIST)) {
            // If this array is the first condition or if its length is equal to the existing length,
            // then it is not necessary to expand the existing data
            if (st) {
                // Overall mask and result arrays have already been established
                // Copy the boolean array
                for (int i=arglen-1; i >= 0; i--) stackpos[st*arglen+i] = stackpos[st+i];
                last_mask_ptr = stackpos + st*arglen;
                // Expand the value data. We can't use memcpy here because if arglen=2 and the result
                // is a point then the point at the lowest index will overlap its destination.
                for (int i=arglen*GET_STEP(result_type)-1; i >= 0; i--) stackpos[arglen+i] = stackpos[1+i%result_step];
                result_ptr = stackpos+arglen;
                st = st*arglen;
            } else if (last_mask_ptr) {
                // First condition has been evaluated, but the first condition is a scalar
                // and the first value is an array.
                // Copy the value
                for (int i=arglen-1; i >= 0; i--) stackpos[i+arglen] = stackpos[i+1];
                lastlen = arglen;
            }
            // Expand the (initial) mask
            for (int i=1; i < arglen; i++) stackpos[i] = stackpos[0];
            result_type |= TYPE_LIST;
            result_len = arglen;
        }
        // Shrink the existing arrays if necessary
        if (argtype & TYPE_LIST) result_len = (arglen/step < result_len ? arglen/step : result_len);
        if ((argtype & TYPE_MASK) == TYPE_BOOLEAN) {
            // Condition was found. Thus, the previous expression was a double and lastlen will be zero.
            last_mask_ptr = stackpos+st;
            lastlen = arglen;
            //printf("    last_mask_ptr %p (%ld)\n", last_mask_ptr, last_mask_ptr - stackpos);
        } else {
            // For the first value, set result_type
            if (!skipped) {
                if ((st) && ((result_type & TYPE_MASK) != (argtype & TYPE_MASK))) FAIL("ERROR: all branches of a conditional must have the same type (%p)\n", fs);
                result_type |= argtype;
                result_step = GET_STEP(result_type);
            }
            value_ptr = stackpos+st+lastlen;
            if (!(argtype & TYPE_LIST)) {
                //printf("casting, value at %p, result_len %d, step %d\n", value_ptr, result_len, step);
                for (int i=1; i < result_len; i++) memcpy(value_ptr+i*step, value_ptr, step*sizeof(double));
            }
            //printf("    After casting: "); print_object((argtype & TYPE_MASK) | TYPE_LIST | ((step*result_len)<<8), value_ptr); printf(", at %p, %d, %d\n", value_ptr, result_len, step);
            if (!st) {
                // If this is the first value, then the combined mask list has not
                // yet been created. We must create the combined mask list and copy
                // the first mask to a later position, after the first value.

                // v- st
                // v- last_mask_ptr
                // | mask_1 | value_1 |

                result_ptr = stackpos+lastlen;
                value_ptr = stackpos+lastlen;
                st = lastlen + step*lastlen;
                last_mask_ptr = stackpos+st;
                memcpy(stackpos+st, stackpos, result_len*sizeof(double));
                for (int i=0; i < result_len; i++) stackpos[i] = 0;
                // | mask | result_1 | mask_1
                //        ^- result_ptr
                //        ^- value_ptr
                //                   ^- st
                //                   ^- last_mask_ptr
            }
            //printf("    After casting: "); print_object((argtype & TYPE_MASK) | TYPE_LIST | ((step*result_len)<<8), value_ptr); printf(", at %p, %d, %d\n", value_ptr, result_len, step);
            if ((lasttype & TYPE_MASK) != TYPE_BOOLEAN) {
                // Default expression
                last_mask_ptr = stackpos+st+lastlen+result_len*step;
                for (int i=0; i < result_len; i++) last_mask_ptr[i] = (stackpos[i] != 0 ? 0 : 1);
                lastlen = result_len;
                //printf("    For default expression, using "); print_object(result_len<<8, last_mask_ptr); printf("\n");
            }
            if (last_mask_ptr) {
                //printf("    last_mask_ptr (len %d) ", lastlen); print_object(lastlen<<8, last_mask_ptr); printf("\n");
                // list of booleans
                if (value_ptr == result_ptr) {
                    for (int i=0; i < result_len; i++) {
                        if (!(last_mask_ptr[i%lastlen])) {
                            for (int j=0; j < step; j++) result_ptr[i*step+j] = NAN;
                        }
                    }
                } else {
                    for (int i=0; i < result_len; i++) {
                        if (last_mask_ptr[i%lastlen] && !stackpos[i]) {
                            //printf("    copying for index %d\n", i);
                            memcpy(result_ptr+i*step, value_ptr+i*step, step*sizeof(double));
                        }
                    }
                }
                // Update the combined mask
                for (int i=0; i < result_len; i++) stackpos[i] += last_mask_ptr[i%lastlen];
            }
            // Short-circuited evaluation
            if ((lastlen == 1) && (last_mask_ptr[0] == 1)) break;
            lastlen = 0;
            // Break if a default expression was found
            if ((lasttype & TYPE_MASK) != TYPE_BOOLEAN) break;
        }
        //printf("combined mask "); print_object(result_len<<8, stackpos);
        //if (result_ptr) {
        //    printf("\nresult "); print_object(((result_len*step)<<8) | (result_type & 0xff), result_ptr);
        //}
        //printf("\n");
        lasttype = argtype;
        arg = arg->next_arg;
    }
    result_len = result_len * GET_STEP(result_type);
    //printf("\nresult "); print_object((result_len<<8) | (result_type & 0xff), result_ptr); printf("\n");
    for (int i=0; i < result_len; i++) stackpos[i] = result_ptr[i]*SIGN_BIT(fs);
    //exit(EXIT_FAILURE);
    return (result_len<<8) | (result_type & 0xff);
}

int compare_doubles(const void *a, const void *b) {
    double arg1 = *(const double*)a;
    double arg2 = *(const double*)b;
    return (arg1 > arg2) - (arg1 < arg2);
}

#ifdef QSORT_MACOS
int compare_doubles_indirect(void *data, const void *a, const void *b) {
#else
int compare_doubles_indirect(const void *a, const void *b, void *data) {
#endif
    double *data_d = (double*)data;
    double arg1 = data_d[(int)(*(double*)a)];
    double arg2 = data_d[(int)(*(double*)b)];
    return (arg1 > arg2) - (arg1 < arg2);
}

uint32_t func_sort(void *f, double *stackpos) {
    // [7,3,5,8,2,1,5,6,8] sorted by [1,7,4,2,6,3]
    // Indices are [0, 3, 5, 2, 4, 1]
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t argtype, idxlen, arglen, result_len, argsize;
    if (arg->next_arg) {
        argtype = arg->next_arg->oper(arg->next_arg, stackpos);
        idxlen = argtype>>8;
        for (int i=0; i < idxlen; i++) {
            stackpos[idxlen+i] = stackpos[i];
            stackpos[i] = i;
        }
#ifdef QSORT_MACOS
        qsort_r(stackpos, idxlen, sizeof(double), stackpos+idxlen, compare_doubles_indirect);
#else
        qsort_r(stackpos, idxlen, sizeof(double), compare_doubles_indirect, stackpos+idxlen);
#endif
        argtype = arg->oper(arg, stackpos+idxlen);
        argsize = argtype>>8;
        uint8_t step = GET_STEP(argtype);
        arglen = argsize/step;
        int j=0;
        for (int i=0; i < idxlen; i++) {
            if (stackpos[i] < arglen) {
                memcpy(stackpos+idxlen+argsize+j, stackpos+idxlen+step*((int)stackpos[i]), step*sizeof(double));
                j += step;
            }
        }
        for (int i=0; i < j; i++) stackpos[i] = stackpos[idxlen+argsize+i]*SIGN_BIT(fs);
        return (j<<8) | (argtype & 0xff);
    } else {
        // Single-argument sort
        argtype = arg->oper(arg, stackpos);
        if (!IS_TYPE(argtype, TYPE_DOUBLE)) FAIL("ERROR: can only sort list of doubles\n");
        qsort(stackpos, argtype>>8, sizeof(double), compare_doubles);
        for (int i=0; i < (argtype>>8); i++) stackpos[i] *= SIGN_BIT(fs);
        return argtype;
    }
}

#ifdef QSORT_MACOS
int compare_generic_indirect(void *data_v, const void *a, const void *b) {
#else
int compare_generic_indirect(const void *a, const void *b, void *data_v) {
#endif
    void **data_s = (void**)(data_v);
    uint8_t step = *((uint8_t*)(data_s[0]));
    double *data = (double*)(data_s[1]);
    double *arg1 = data + (step*(int)(*((double*)a)));
    double *arg2 = data + (step*(int)(*((double*)b)));
    for (int i=0; i < step; i++) {
        if (arg1[i] > arg2[i]) return 1;
        if (arg1[i] < arg2[i]) return -1;
    }
    return 0;
}

uint32_t func_unique(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t argtype = arg->oper(arg, stackpos);
    uint32_t arglen = argtype>>8;
    uint8_t step = GET_STEP(argtype);
    uint32_t n = arglen / step;
    for (int i=0; i < n; i++) stackpos[arglen+i] = i;
    void *data[2] = {&step, stackpos};
#ifdef QSORT_MACOS
    qsort_r(stackpos+arglen, n, sizeof(double), data, compare_generic_indirect);
#else
    qsort_r(stackpos+arglen, n, sizeof(double), compare_generic_indirect, data);
#endif
    int j=arglen+1;
    for (int i=1; i < n; i++) {
        //print_object(256 | TYPE_POINT, stackpos+step*(int)stackpos[j-1]); printf(", "); print_object(256 | TYPE_POINT, stackpos+step*(int)stackpos[arglen+i]); printf("\n");
        if (memcmp(stackpos+step*(int)stackpos[j-1], stackpos+step*(int)stackpos[arglen+i], step*sizeof(double))) {
            stackpos[j] = stackpos[arglen+i];
            j++;
        }
    }
    j -= arglen;
    qsort(stackpos+arglen, j, sizeof(double), compare_doubles);
    for (int i=0; i < j; i++) {
        if (i < stackpos[arglen+i]) memcpy(stackpos+i*step, stackpos+step*(int)stackpos[arglen+i], step*sizeof(double));
    }
    return ((j*step)<<8) | (argtype & 0xff);
}

uint32_t func_join(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    
    uint32_t st=0;
    uint32_t type=-1, argtype;
    while (arg) {
        argtype = arg->oper(arg, stackpos+st);
        if ((type != -1) && ((type & TYPE_MASK) != (argtype & TYPE_MASK))) FAIL("ERROR: lists must have same type to join\n");
        type = argtype;
        st += argtype>>8;
        arg = arg->next_arg;
    }
    return (st<<8) | (type & 0xff) | TYPE_LIST;
}

uint32_t func_length(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    
    uint32_t argtype = arg->oper(arg, stackpos);
    stackpos[0] = (argtype>>8)/GET_STEP(argtype);
    return 1<<8;
}

double integrate_func(double x, void *params) {
    function *expr = (function*)(((void**)params)[0]);
    double *stackpos = (double*)(((void**)params)[1]);
    variable *var = (variable*)(((void**)params)[2]);
    var->pointer = &x;
    var->type = 1<<8;
    expr->oper(expr, stackpos);
    var->pointer = 0;
    return stackpos[0];
}

uint32_t func_integrate_gsl(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *var = fs->first_arg;
    function *lb = var->next_arg;
    function *ub = lb->next_arg;
    function *expr = ub->next_arg;

    lb->oper(lb, stackpos);
    double lbv = stackpos[0];
    ub->oper(ub, stackpos);
    double ubv = stackpos[0];
    integration_result *res = (integration_result*)(fs->value);
    if (lbv == ubv) {
        // Null case is not cached
        stackpos[0] = 0;
        return 1<<8;
    }
    uint8_t negate = 0;
    if (lbv > ubv) {
        negate = 1;
        double t = lbv;
        lbv = ubv;
        ubv = t;
    }
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    void *params[3] = {expr, stackpos, var->value};
    F.function = &integrate_func;
    F.params = params;
    double result, error;
    gsl_set_error_handler_off();
    int status;
    if (res && (res->func == fs) && (ubv-lbv >= res->ubv-res->lbv)) {
        if (lbv == res->lbv) {
            if (res->result == res->result) status = gsl_integration_qags(&F, res->ubv, ubv, 0, 1e-7, 1000, w, &result, &error);
            result += res->result;
        } else if (ubv == res->ubv) {
            if (res->result == res->result) status = gsl_integration_qags(&F, lbv, res->lbv, 0, 1e-7, 1000, w, &result, &error);
            result += res->result;
        } else status = gsl_integration_qags(&F, lbv, ubv, 0, 1e-7, 1000, w, &result, &error);
    } else status = gsl_integration_qags(&F, lbv, ubv, 0, 1e-7, 1000, w, &result, &error);
    stackpos[0] = (negate ? -result : result);
    if (status) {
        printf("ERROR: error %d while integrating %p from %f to %f, value %f: %s\n", status, expr, lbv, ubv, result, gsl_strerror(status));
        if (res) {
            printf("cached %f from %f to %f\n", res->result, res->lbv, res->ubv);
        }
        //stackpos[0] = NAN;
        //result = NAN;
    }
    /*else if (status) {
        FAIL("ERROR: error %d while integrating from %f to %f, value %f: %s\n", status, lbv, result, ubv, gsl_strerror(status));
    }*/
    if (res) {
        res->lbv = lbv;
        res->ubv = ubv;
        res->result = result;
        res->func = fs;
        res->negate = negate;
    }
    gsl_integration_workspace_free(w);
    //printf("integrate from %f to %f, result %f\n", lbv, ubv, stackpos[0]);
    return 1<<8;
}

double integrate_rec(function *func, double *stackpos, double *var, double lbv, double ubv, double *buf) {
    double result = 0;
    double xpos = INTEGRATE_DX;
    if (lbv > ubv) return -integrate_rec(func, stackpos, var, ubv, lbv, buf);
    *var = lbv;
    func->oper(func, stackpos);
    double yp = stackpos[0];
    if (buf) {
        buf[0] = *var;
        buf[1] = yp;
    }
    uint32_t i=2;
    while (*var < ubv) {
        func->oper(func, stackpos);
        result += INTEGRATE_DX*(stackpos[0] + yp)/2;
        yp = stackpos[0];
        if (buf) {
            buf[i++] = *var;
            buf[i++] = yp;
        }
        xpos += INTEGRATE_DX;
        *var = lbv + xpos;
    }
    // final step
    *var = ubv;
    func->oper(func, stackpos);
    if (buf) {
        buf[i++] = *var;
        buf[i++] = stackpos[0];
    }
    result += (ubv - (lbv + (xpos - INTEGRATE_DX)))*(stackpos[0] + yp)/2;
    return result;
}

uint32_t func_convert_polar(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *radius = fs->first_arg;
    function *angle = radius->next_arg;

    uint32_t atype = angle->oper(angle, stackpos);
    double sine = sin(stackpos[0]);
    double cosine = cos(stackpos[0]);
    uint32_t rtype = radius->oper(radius, stackpos);
    uint32_t len = rtype>>8;

    for (int i=len; i >= 0; i--) {
        stackpos[2*i+1] = stackpos[i]*sine;
        stackpos[2*i] = stackpos[i]*cosine;
    }
    return (len<<9) | (rtype & TYPE_LIST) | TYPE_POINT;
}

uint32_t func_onkeypress(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *key = fs->first_arg;
    function *action = key->next_arg;
    key->oper(key, stackpos);
    uint32_t code = (uint32_t)(stackpos[0]);
    printf("Assigning %p to keypress event %d\n", action, code);
    *((uint64_t*)(stackpos+1)) = (uint64_t)(action);
    return 1<<9;
}

uint32_t func_color(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *target = fs->first_arg;
    return target->oper(target, stackpos);
}

#define N_OPERATORS 52
const oper_data oper_list[N_OPERATORS] = {
    {NULL, "unknown_function", NULL},
    {func_value, "func_value", interval_value},

    {func_div, "func_div", interval_div},
    {func_floor, "func_floor", NULL},
    {func_round, "func_round", NULL},
    {func_mod, "func_mod", NULL},
    {func_max, "func_max", NULL},
    {func_min, "func_min", NULL},
    {func_sign, "func_sign", NULL},
    {func_sine, "func_sine", interval_sine},
    {func_cosine, "func_cosine", interval_cosine},
    {func_tangent, "func_tangent", interval_tangent},
    {func_arctan, "func_arctan", interval_arctan},
    {func_arcsin, "func_arcsin", NULL},
    {func_conjugate, "func_conjugate", interval_conjugate},
    {func_log, "func_log", NULL},
    {func_exp, "func_exp", interval_exp},
    {func_factorial, "func_factorial", interval_factorial},
    {func_abs, "func_abs", interval_abs},
    {func_add, "func_add", interval_add},
    {func_multiply, "func_multiply", interval_multiply},
    {func_exponentiate, "func_exponentiate", interval_exponentiate},
    {func_user_defined, "func_user_defined", interval_user_defined},

    {func_list, "func_list", interval_list},
    {func_index, "func_index", interval_index},
    {func_point, "func_point", interval_point},
    {func_polygon, "func_polygon", NULL},
    {func_rgb, "func_rgb", NULL},
    {func_ellipsis, "func_ellipsis", NULL},

    {func_for, "func_for", NULL},
    {func_equals, "func_equals", interval_equals},
    {func_greater, "func_greater", interval_greater},
    {func_compare_single, "func_compare_single", interval_compare_single},
    {func_compare, "func_compare", interval_compare},

    {func_extract_x, "func_extract_x", interval_extract_x},
    {func_extract_y, "func_extract_y", interval_extract_y},
    {func_assign, "func_assign", NULL},
    {func_chain_actions, "func_chain_actions", NULL},
    {func_sum, "func_sum", NULL},
    {func_prod, "func_prod", NULL},
    {func_total, "func_total", NULL},
    {func_distance, "func_distance", NULL},
    {func_conditional, "func_conditional", NULL},
    {func_sort, "func_sort", interval_sort},
    {func_unique, "func_unique", NULL},
    {func_join, "func_join", NULL},
    {func_length, "func_length", NULL},
    {func_integrate_gsl, "func_integrate_gsl", interval_integrate_gsl},
    
    {func_det, "func_det", interval_det},

    {func_convert_polar, "func_convert_polar", NULL},
    {func_onkeypress, "func_onkeypress", NULL},
    {func_color, "func_color", NULL}
};

const oper_data *oper_lookup(uint32_t (*ptr)(void*, double*)) {
    for (int i=0; i < N_OPERATORS; i++) {
        if (oper_list[i].oper == ptr) return oper_list+i;
    }
    return oper_list;
}
