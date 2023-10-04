#include "functions.h"
#include "parse.h"
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SIGN_BIT(v) (((v->value_type)&0x80) ? -1 : 1)
#define VALUE(a) ((((a->value_type)&0x40) ? *(((variable*)(a->value))->pointer) : *((double*)(a->value))) * SIGN_BIT(a))
#define VALUE_LIST(a, idx) ((((a->value_type)&0x40) ? ((double*)(((variable*)(a->value))->pointer))[idx] : *((double*)(a->value) + idx)) * SIGN_BIT(a))

#define FAIL(...) {printf(__VA_ARGS__); exit(EXIT_FAILURE);}

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
    arg = arg->next_arg;
    if (arg->oper) {
        val2 = stackpos+(l1?l1:1);
        t2 = arg->oper(arg, val2);
    } else {
        t2 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
        val2 = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->pointer : (double*)(arg->value));
    }
    l1 = t1>>8;
    l2 = t2>>8;
    if (!(t1 & TYPE_LIST) && !(t2 & TYPE_LIST)) {
        stackpos[0] = oper(*val1, *val2)*SIGN_BIT(fs);
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

double fdiv(double n, double d) {
    return n/d;
}

uint32_t func_div(void *f, double *stackpos) {
    return func_general_two_args(f, stackpos, fdiv);
}

uint32_t func_sine(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, sin);
}

uint32_t func_cosine(void *f, double *stackpos) {
    return func_general_one_arg(f, stackpos, cos);
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
                sum += stackpos[0];
            }
        }
        arg = arg->next_arg;
    } while (arg);
    if (!(result_type & TYPE_LIST)) {
        stackpos[0] = sum*SIGN_BIT(fs);
        result_length = 1;
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
                prod *= stackpos[0];
            }
        }
        arg = arg->next_arg;
    } while (arg);
    if (!(result_type & TYPE_LIST)) {
        stackpos[0] = prod*SIGN_BIT(fs);
        result_length = 1;
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
    do {
        if (arg->oper) {
            type = arg->oper(arg, stackpos+st);
            target_arg[varnum].pointer = stackpos+st;
            varnum++;
            st += (type>>8);
        } else {
            type = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
            len = type>>8;
            for (int i=0; i < len; i++) stackpos[st+i] = VALUE_LIST(arg, i);
            target_arg[varnum].pointer = stackpos+st;
            target_arg[varnum].type = type;
            varnum++;
            st += len;
        }
        arg = arg->next_arg;
    } while (arg);
    type = target->oper(target, stackpos+st);
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
    do {
        if (arg->oper) {
            argtype = arg->oper(arg, stackpos+st);
            if (argtype & TYPE_LIST) FAIL("ERROR: Cannot store a list inside a list\n");
            if ((type != 0xff) && ((argtype&0xff) != type)) FAIL("ERROR: All list elements must have the same type\n");
            type = argtype & 0xff;
            st += argtype>>8;
        } else {
            argtype = ((arg->value_type)&0x40 ? ((variable*)(arg->value))->type : arg->value_type);
            //printf("First argument: %08x, %f\n", argtype, VALUE(arg));
            if (argtype & TYPE_LIST) FAIL("ERROR: Cannot store a list inside a list\n");
            if ((type != 0xff) && ((argtype&0xff) != type)) FAIL("ERROR: All list elements must have the same type\n");
            type = argtype & 0xff;
            arglen = argtype>>8;
            for (int i=0; i < arglen; i++) stackpos[st+i] = VALUE_LIST(arg, i);
            st += arglen;
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
        int32_t v2 = step*(val1[0]-1);
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
        return (l1<<9) | TYPE_POINT;
    } else if (!(t1 & TYPE_LIST) && (t2 & TYPE_LIST)) {
        double v1 = *val1;
        for (int i=l2-1; i >= 0; i--) {
            stackpos[2*i] = v1;
            stackpos[2*i+1] = stackpos[i+1]*SIGN_BIT(fs);
        }
        return (l2<<9) | TYPE_POINT;
    } else if (l1 > l2) {
        memcpy(stackpos+l1+l2, stackpos, l2);
        val1 = stackpos+l1+l2;
        for (int i=0; i < l2; i++) {
            stackpos[2*i] = val1[i]*SIGN_BIT(fs);
            stackpos[2*i+1] = val2[i]*SIGN_BIT(fs);
        }
        return (l2<<9) | TYPE_POINT;
    } else {
        memcpy(stackpos+l1+l2, stackpos, l1);
        val1 = stackpos+l1+l2;
        for (int i=0; i < l1; i++) {
            stackpos[2*i] = val1[i]*SIGN_BIT(fs);
            stackpos[2*i+1] = val2[i]*SIGN_BIT(fs);
        }
        return (l1<<9) | TYPE_POINT;
    }
    
}

