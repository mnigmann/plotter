#include "functions.h"
#include "parse.h"
#include "config.h"
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SIGN_BIT(v) (((v->value_type)&0x80) ? -1 : 1)
#define VALUE(a) ((((a->value_type)&0x40) ? *(((variable*)(a->value))->pointer) : *((double*)(a->value))) * SIGN_BIT(a))
#define VALUE_LIST(a, idx) ((((a->value_type)&0x40) ? ((double*)(((variable*)(a->value))->pointer))[idx] : *((double*)(a->value) + idx)) * SIGN_BIT(a))

#define FAIL(...) {printf(__VA_ARGS__); exit(EXIT_FAILURE);}

const double lanczos_table[9] = {
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};


uint32_t func_general_one_arg(void *f, double *stackpos, double (*oper)(double)) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t len;
    if (arg->oper) {
        len = arg->oper(arg, stackpos)>>8;
        if (len) {
            for (int i=0; i < len; i++) stackpos[i] = oper(stackpos[i])*SIGN_BIT(fs);
        } else stackpos[0] = oper(stackpos[0])*SIGN_BIT(fs);
        return len<<8;
    } else {
        len = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type)>>8;
        if (len) {
            for (int i=0; i < len; i++) stackpos[i] = oper(VALUE_LIST(arg, i)) * SIGN_BIT(fs);
        } else stackpos[0] = oper(VALUE(arg)) * SIGN_BIT(fs);
        return arg->value_type;
    }
}

uint32_t func_general_two_args(void *f, double *stackpos, double (*oper)(double, double)) {
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
    l1 = t1>>8;
    arg = arg->next_arg;
    if (arg->oper) {
        val2 = stackpos+(l1?l1:1);
        t2 = arg->oper(arg, val2);
    } else {
        t2 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
        val2 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->pointer : (double*)(arg->value));
    }
    l2 = t2>>8;
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
    uint32_t type = ((fs->value_type)&0x40 ? ((variable*)(fs->value))->type : fs->value_type);
    double *ptr = ((fs->value_type)&0x40 ? ((variable*)(fs->value))->pointer : (double*)(fs->value));
    uint32_t len = type>>8;
    for (int i=0; i < len; i++) stackpos[i] = ptr[i]*SIGN_BIT(fs);
    return type;
}

double fdiv(double n, double d) {
    return n/d;
}

double mmod(double n, double d) {
    double i = fmod(n, d);
    if (n && ((n>0) ^ (d>0))) i += d;
    return i;
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

uint32_t func_div(void *f, double *stackpos) {
    return func_general_two_args(f, stackpos, fdiv);
}

uint32_t func_mod(void *f, double *stackpos) {
    return func_general_two_args(f, stackpos, mmod);
}

uint32_t func_floor(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, floor);
}

uint32_t func_sine(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, sin);
}

uint32_t func_cosine(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, cos);
}

uint32_t func_factorial(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, mfac);
}

// The arctan function can have either one or two arguments
uint32_t func_arctan(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    if (arg->next_arg) return func_general_two_args(f, stackpos, atan2);
    else return func_general_one_arg(f, stackpos, atan);
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
    do {
        if (arg->oper) {
            type = arg->oper(arg, stackpos+result_length);
            len = type>>8;
        } else {
            type = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
            len = type>>8;
            for (i=0; i < len; i++) stackpos[i+result_length] = VALUE_LIST(arg, i);
        }
        if (type & TYPE_LIST) {
            if (result_type & TYPE_LIST) {
                // The result has already been made a list, so we multiply and pick the minimum length
                for (i=0; (i < result_length) && (i < len); i++) stackpos[i] += stackpos[result_length+i];
                // i will be the minimum of result_length and len
                result_length = i;
                result_type = TYPE_LIST;
            } else {
                // The result is currently a scalar and must be updated
                for (i=0; i < len; i++) stackpos[i] += sum;
                result_length = len;
                result_type = TYPE_LIST;
            }
        } else {
            if (result_type & TYPE_LIST) {
                // The result is already a list but the current term is a scalar
                for (i=0; i < result_length; i++) stackpos[i] += stackpos[result_length];
            } else {
                // The result is a scalar and the term is a scalar
                sum += stackpos[result_length];
            }
        }
        arg = arg->next_arg;
    } while (arg);
    if (!(result_type & TYPE_LIST)) {
        stackpos[0] = sum*SIGN_BIT(fs);
        result_length = 1;
    } else {
        for (i=0; i < result_length; i++) stackpos[i] *= SIGN_BIT(fs);
    }
    return (result_length<<8) | result_type;
}

uint32_t func_multiply(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    double prod = 1;
    uint32_t type;
    uint32_t len;
    uint32_t result_type = 0;
    uint32_t result_length = 0;
    int i=0;
    do {
        if (arg->oper) {
            type = arg->oper(arg, stackpos+result_length);
            len = type>>8;
        } else {
            type = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
            len = type>>8;
            for (i=0; i < len; i++) stackpos[i+result_length] = VALUE_LIST(arg, i);
        }
        if (type & TYPE_LIST) {
            if (result_type & TYPE_LIST) {
                // The result has already been made a list, so we multiply and pick the minimum length
                for (i=0; (i < result_length) && (i < len); i++) stackpos[i] *= stackpos[result_length+i];
                // i will be the minimum of result_length and len
                result_length = i;
                result_type = TYPE_LIST;
            } else {
                // The result is currently a scalar and must be updated
                for (i=0; i < len; i++) stackpos[i] *= prod;
                result_length = len;
                result_type = TYPE_LIST;
            }
        } else {
            if (result_type & TYPE_LIST) {
                // The result is already a list but the current term is a scalar
                for (i=0; i < result_length; i++) stackpos[i] *= stackpos[result_length];
            } else {
                // The result is a scalar and the term is a scalar
                prod *= stackpos[result_length];
            }
        }
        arg = arg->next_arg;
    } while (arg);
    if (!(result_type & TYPE_LIST)) {
        stackpos[0] = prod*SIGN_BIT(fs);
        result_length = 1;
    } else {
        for (i=0; i < result_length; i++) stackpos[i] *= SIGN_BIT(fs);
    }
    return (result_length<<8) | result_type;
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
        //printf("argument is %08x:", type); print_object(type, target_arg[varnum].pointer);
        varnum++;
        arg = arg->next_arg;
    } while (arg);
    type = target->oper(target, stackpos+st);
    //printf("return type is %08x\n", type);
    len = type>>8;
    // Shrink the stack by removing unused values
    for (int i=0; i < len; i++) stackpos[i] = stackpos[st+i]*SIGN_BIT(fs);
    return type;
}

uint32_t func_list(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint8_t type = 0xff;
    uint32_t argtype;
    uint32_t arglen;
    uint32_t st=0;
    uint8_t ellipsis = 0;
    do {
        if (arg->oper) {
            argtype = arg->oper(arg, stackpos+st);
            if ((argtype == TYPE_ELLIPSIS) && ((type & TYPE_MASK) == TYPE_DOUBLE) && (arg->next_arg)) {
                ellipsis = 2;
            } else {
                if (argtype & TYPE_LIST) FAIL("ERROR: Cannot store a list inside a list\n");
                if ((type != 0xff) && ((argtype&0xff) != type)) FAIL("ERROR: All list elements must have the same type\n");
                st += argtype>>8;
                type = argtype & 0xff;
            }
        } else {
            argtype = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
            if (argtype & TYPE_LIST) FAIL("ERROR: Cannot store a list inside a list\n");
            if ((type != 0xff) && ((argtype&0xff) != type)) FAIL("ERROR: All list elements must have the same type\n");
            type = argtype & 0xff;
            arglen = argtype>>8;
            for (int i=0; i < arglen; i++) stackpos[st+i] = VALUE_LIST(arg, i);
            st += arglen;
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
        }
        arg = arg->next_arg;
    } while (arg);
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
    if (arg->oper) {
        t1 = arg->oper(arg, stackpos);
        l1 = t1 >> 8;
    } else {
        t1 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
        l1 = t1 >> 8;
        for (int i=0; i < l1; i++) stackpos[i] = VALUE_LIST(arg, i);
    }
    val1 = stackpos;
    if (!(t1 & TYPE_LIST)) FAIL("ERROR: Cannot index a non-list\n");
    
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

    uint8_t step=1;
    if ((t1 & TYPE_MASK) == TYPE_POINT) step=2;
    if ((t1 & TYPE_MASK) == TYPE_COLOR) step=3;
    if ((t2 & TYPE_LIST) && ((t2 & TYPE_MASK) == TYPE_DOUBLE)) {
        // List of indices
        int32_t v2;
        for (int i=0; i < l2; i++) {
            v2 = step*((int32_t)val2[i]-1);
            for (uint8_t j=0; j < step; j++) stackpos[l1+i*step+j] = (((v2 >= l1) || (v2 < 0)) ? NAN : val1[v2+j]);
        }
        for (int i=0; i < l2*step; i++) stackpos[i] = stackpos[l1+i];
        return ((step*l2)<<8) | (t1 & 0xff);
    } else if ((t2 & TYPE_LIST) && ((t2 & TYPE_MASK) == TYPE_BOOLEAN)) {
        // List of booleans
        uint32_t lmin = (l1 < l2 ? l1 : l2);
        uint32_t lr = 0;
        for (int i=0; i < lmin; i++) {
            if (val2[i] == 0) continue;
            for (uint8_t j=0; j < step; j++) stackpos[lr+j] = val1[i*step+j];
            lr+=step;
        }
        return (lr << 8) | (t1 & 0xff);
    } else if ((t2 & TYPE_MASK) == TYPE_DOUBLE) {
        int32_t v2 = step*(val2[0]-1);
        for (uint8_t j=0; j < step; j++) stackpos[j] = (((v2 >= l1) || (v2 < 0)) ? NAN : val1[v2+j]);
        return (step << 8) | (t2 & 0xff & ~TYPE_LIST);
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
        for (int j=0; j < MAX_POLYGON_SIZE; j++) {
            stackpos[i*2*MAX_POLYGON_SIZE + 2*j] = temp[j*first_len + 2*i];
            stackpos[i*2*MAX_POLYGON_SIZE + 2*j + 1] = temp[j*first_len + 2*i + 1];
        }
    }
    free(temp);
    return ((min_len*MAX_POLYGON_SIZE)<<8) | TYPE_LIST | TYPE_POLYGON;
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
    for (i=input_len; i < st; i++) stackpos[i-input_len] = stackpos[i];
    
    return ((st - input_len)<<8) | (0xff&argtype) | TYPE_LIST;
}

uint32_t func_equals(void *f, double *stackpos) {
    return TYPE_BOOLEAN;
}

uint32_t func_extract_x(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    //printf("extracting x, arg %p, oper %p (%p), value %p, value_type %08x\n", arg, arg->oper, (func_value), arg->value, arg->value_type);
    
    uint32_t argtype, arglen;
    if (arg->oper) {
        argtype = arg->oper(arg, stackpos);
        arglen = argtype>>8;
    } else {
        argtype = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
        arglen = argtype>>8;
        for (int j=0; j < arglen; j++) stackpos[j] = VALUE_LIST(arg, j);
    }
    //printf("extracting from type %08x\n", argtype);
    if (!((argtype & TYPE_MASK) == TYPE_POINT)) FAIL("ERROR: cannot access x-coordinate of type %08x\n", argtype);
    for (int i=0; i < arglen; i += 2) {
        stackpos[i>>1] = stackpos[i];
    }
    return ((argtype>>1) & 0xffffff00) | ((argtype&0xff) & ~TYPE_MASK);
}

uint32_t func_extract_y(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    
    uint32_t argtype, arglen;
    if (arg->oper) {
        argtype = arg->oper(arg, stackpos);
        arglen = argtype>>8;
    } else {
        argtype = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
        arglen = argtype>>8;
        for (int j=0; j < arglen; j++) stackpos[j] = VALUE_LIST(arg, j);
    }
    if (!((argtype & TYPE_MASK) == TYPE_POINT)) FAIL("ERROR: cannot access y-coordinate of type %08x\n", argtype);
    for (int i=0; i < arglen; i += 2) {
        stackpos[i>>1] = stackpos[i+1];
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

    // Call the next action in the list. Comma-separated sequences are automatically
    // parsed into chained blocks
    if (fs->next_arg) fs->next_arg->oper(fs->next_arg, stackpos);

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
    
    uint32_t result_length = 0;
    uint32_t result_type = 0;
    double value;
    double sum = 0;
    int i;
    double start_val = *val1;
    double end_val = *val2;
    if (!(t_start & TYPE_LIST) && !(t_end & TYPE_LIST)) {
        varptr->pointer = &value;
        varptr->type = 0x100;
        for (value=start_val; value <= end_val; value++) {
            t_expr = expr->oper(expr, stackpos+result_length);
            uint32_t len = t_expr>>8;
            if (t_expr & TYPE_LIST) {
                if (result_type & TYPE_LIST) {
                    // The result has already been made a list, so we multiply and pick the minimum length
                    for (i=0; (i < result_length) && (i < len); i++) stackpos[st+i] += stackpos[result_length+i];
                    // i will be the minimum of result_length and len
                    result_length = i;
                    result_type = TYPE_LIST;
                } else {
                    // The result is currently a scalar and must be updated
                    for (i=0; i < len; i++) stackpos[i] += sum;
                    result_length = len;
                    result_type = TYPE_LIST;
                }
            } else {
                if (result_type & TYPE_LIST) {
                    // The result is already a list but the current term is a scalar
                    for (i=0; i < result_length; i++) stackpos[i] += stackpos[result_length];
                } else {
                    // The result is a scalar and the term is a scalar
                    sum += stackpos[result_length];
                }
            }
        }
    }
    if (!(result_type & TYPE_LIST)) {
        stackpos[0] = sum*SIGN_BIT(fs);
        result_length = 1;
    } else {
        for (int j=0; j < result_length; j++) stackpos[j] = stackpos[j]*SIGN_BIT(fs);
    }
    return (result_length<<8) | result_type;
}

