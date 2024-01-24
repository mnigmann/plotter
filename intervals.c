#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "intervals.h"
#include "parse.h"
#include "functions.h"

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

#define MAX_POLYGON_SIZE 4
const static uint32_t step_table[8] = {1, 2, 3, 1, 1, MAX_POLYGON_SIZE*2, 0, 0};

/*
 * Generic helper functions
 */
void apply_sign(uint8_t sign_bit, double *hstack, double *lstack, uint32_t srcpos, uint32_t len) {
    if (sign_bit & 0x80) {
        double temp;
        for (uint32_t i=0; i < len; i++) {
            temp = hstack[i+srcpos];
            hstack[i] = -lstack[i+srcpos];
            lstack[i] = -temp;
        }
    } else if (srcpos) {
        for (uint32_t i=0; i < len; i++) {
            hstack[i] = hstack[i+srcpos];
            lstack[i] = lstack[i+srcpos];
        }
    }
}

uint32_t interval_general_one_arg(void *f, double *hstackpos, double *lstackpos, void (*oper)(double*, double*)) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t len, type;
    
    type = arg->inter(arg, hstackpos, lstackpos);
    len = type>>8;
    for (int i=0; i < len; i++) oper(hstackpos+i, lstackpos+i);
    apply_sign(fs->value_type, hstackpos, lstackpos, 0, len);
    return (len<<8) | (type & TYPE_LIST);
}

uint32_t interval_general_two_args(void *f, double *hstackpos, double *lstackpos, void (*oper)(double*, double*, double, double)) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t t1, t2, l1, l2;
    t1 = arg->inter(arg, hstackpos, lstackpos);
    l1 = t1>>8;
    arg = arg->next_arg;
    t2 = arg->inter(arg, hstackpos+l1, lstackpos+l1);
    l2 = t2>>8;
    //printf("func_general_two_args: "); print_object(t1, val1); printf(" and "); print_object(t2, val2); printf("\n");
    uint32_t type;
    uint32_t result_len = l1;
    uint32_t result_pos = 0;
    if (!(t1 & TYPE_LIST) || ((t2 & TYPE_LIST) && (l2 < result_len))) result_len = l2;
    if (!(t1 & TYPE_LIST) && (t2 & TYPE_LIST)) {
        double vh = hstackpos[0], vl = lstackpos[0];
        for (int i=0; i < l2; i++) {
            hstackpos[i] = vh;
            lstackpos[i] = vl;
            oper(hstackpos+i, lstackpos+i, hstackpos[i+1], lstackpos[i+1]);
        }
    } else {
        for (int i=0; i < result_len; i++) {
            oper(hstackpos+i, lstackpos+i, hstackpos[l1+i%l2], lstackpos[l1+i%l2]);
        }
    }
    apply_sign(fs->value_type, hstackpos, lstackpos, 0, result_len);
    return (result_len<<8) | TYPE_LIST;
}

uint32_t interval_value(void *f, double *hstackpos, double *lstackpos) {
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
        if (var->flags & VARIABLE_ACTION) FAIL("ERROR: interval calculation not supported for actions\n");
        double *hpos;
        uint32_t len = type>>8;
        if (var->flags & VARIABLE_INTERVAL) hpos = ptr+len;
        else hpos = ptr;
        if (fs->value_type & 0x80) {
            for (int i=0; i < len; i++) {
                hstackpos[i] = -ptr[i];
                lstackpos[i] = -hpos[i];
            }
        } else {
            for (int i=0; i < len; i++) {
                hstackpos[i] = hpos[i];
                lstackpos[i] = ptr[i];
            }
        }
        //printf("variable %s (%p), flags %02x has interval [%f, %f], len %d, result [%f, %f]\n", var->name, var, var->flags, ptr[0], hpos[0], len, lstackpos[0], hstackpos[0]);
        return type;
    } else {
        type = fs->value_type;
        ptr = (double*)(fs->value);
        if (!ptr) FAIL("ERROR: pointer in function block %p is null\n", fs);
        uint32_t len = type>>8;
        for (int i=0; i < len; i++) {
            hstackpos[i] = ptr[i]*SIGN_BIT(fs);
            lstackpos[i] = hstackpos[i];
        }
        return type;
    }
}

/*
 * Basic arithmetic functions
 */

uint8_t interval_add_in_place(double *hacc, double *lacc, uint32_t *result_type, uint32_t *result_length, uint32_t type) {
    uint32_t len = type>>8;
    uint32_t old_len = *result_length;
    if ((*result_type & TYPE_MASK) != (type & TYPE_MASK)) {
        // Error
        return 1;
    }
    //printf("Adding (%08x): ", type); print_object(type, val); printf(" to "); print_object(*result_type, acc); printf("\n");
    if (!(*result_type & TYPE_LIST) && (type & TYPE_LIST)) {
        // To add a list to a scalar, we swap the operands and then copy down.
        uint8_t step = GET_STEP(*result_type);
        for (int i=0; i < len; i++) {
            hacc[i+old_len] += hacc[i%step];
            lacc[i+old_len] += lacc[i%step];
        }
        *result_length = len;
        *result_type = type;
        for (int i=0; i < len; i++) {
            hacc[i] = hacc[i+old_len];
            lacc[i] = lacc[i+old_len];
        }
        //printf("result "); print_object(*result_type, acc); printf("\n");
    } else {
        // Otherwise, we broadcast and add
        if ((type & TYPE_LIST)  && (len < *result_length)) *result_length = len;
        for (int i=0; i < *result_length; i++) {
            hacc[i] += hacc[i%len + old_len];
            lacc[i] += lacc[i%len + old_len];
        }
    }
    return 0;
}

uint32_t interval_add(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    double sum = 0;
    uint32_t type;
    uint32_t len;
    uint32_t result_type = 0;
    uint32_t result_length = 0;
    int i=0;
    result_type = arg->inter(arg, hstackpos, lstackpos);
    result_length = result_type>>8;
    arg = arg->next_arg;
    double *pos;
    //printf("func_add\n");
    while (arg) {
        type = arg->inter(arg, hstackpos+result_length, lstackpos+result_length);
        //printf("adding %08x: ", type); print_object(type, stackpos+result_length); printf("\n");
        // Keep track of the position in case result_length is changed.
        if (interval_add_in_place(hstackpos, lstackpos, &result_type, &result_length, type)) {
            printf("ERROR: only objects with the same type may be added (function block %p)\n", fs);
            printf("\n    expr:   "); print_object(type, pos);
            FAIL("\n");
        }
        arg = arg->next_arg;
    }
    apply_sign(fs->value_type, hstackpos, lstackpos, 0, result_length);
    return (result_length << 8) | (result_type & 0xff);
}

void single_multiply(double *hstackpos, double *lstackpos, double h, double l) {
    //printf("single_multiply (%f, %f) * (%f, %f)\n", hstackpos[0], lstackpos[0], h, l);
    double ac = lstackpos[0]*l;
    double ad = lstackpos[0]*h;
    double bc = hstackpos[0]*l;
    double bd = hstackpos[0]*h;
    lstackpos[0] = ac;
    hstackpos[0] = ac;
    if (ad > hstackpos[0]) hstackpos[0] = ad;
    if (ad < lstackpos[0]) lstackpos[0] = ad;
    if (bc > hstackpos[0]) hstackpos[0] = bc;
    if (bc < lstackpos[0]) lstackpos[0] = bc;
    if (bd > hstackpos[0]) hstackpos[0] = bd;
    if (bd < lstackpos[0]) lstackpos[0] = bd;
}

void interval_multiply_in_place(double *hstackpos, double *lstackpos, uint32_t *result_type, uint32_t *result_length, uint32_t type) {
    uint32_t len = type>>8;
    uint8_t step = GET_STEP(*result_type);
    uint32_t old_len = *result_length*step;
    if ((!(*result_type & TYPE_LIST) || (len/GET_STEP(type) < *result_length)) && (type & TYPE_LIST)) *result_length = len/GET_STEP(type);
    if (IS_TYPE(*result_type, TYPE_POINT) && IS_TYPE(type, TYPE_POINT)) {
        //printf("multiplying complex\n");
        if (!IS_TYPE(*result_type, TYPE_LIST) && IS_TYPE(type, TYPE_LIST)) {
            // Result is scalar, must convert
            double ah = hstackpos[0], al = lstackpos[0], bh = hstackpos[1], bl = lstackpos[1];
            for (int i=0; i < len; i+=2) {
                hstackpos[i] = ah; lstackpos[i] = al; hstackpos[i+1] = bh; lstackpos[i+1] = bl;
                hstackpos[i] = hstackpos[i+2]; lstackpos[i] = lstackpos[i+2];       // c
                hstackpos[i+1] = hstackpos[i+3]; lstackpos[i+1] = lstackpos[i+3];   // d
                // a*c
                single_multiply(hstackpos+i, lstackpos+i, ah, al);
                // b*d
                single_multiply(hstackpos+i+1, lstackpos+i+1, bh, bl);
                // a*c - d*b
                hstackpos[i] -= lstackpos[i+1];
                lstackpos[i] -= hstackpos[i+1];
                // a*d
                single_multiply(hstackpos+i+3, lstackpos+i+3, ah, al);
                // b*c
                single_multiply(hstackpos+i+2, lstackpos+i+2, bh, bl);
                // a*d + b*c
                hstackpos[i+1] = hstackpos[i+2] + hstackpos[i+3];
                lstackpos[i+1] = lstackpos[i+1] + lstackpos[i+3];
            }
        } else {
            len = len / 2;
            double ah, al, bh, bl;
            for (int i=0; i<*result_length; i++) {
                ah = hstackpos[2*i], al = lstackpos[2*i], bh = hstackpos[2*i+1], bl = lstackpos[2*i+1];
                // a*c
                single_multiply(hstackpos+2*i, lstackpos+2*i, hstackpos[2*i+old_len], lstackpos[2*i+old_len]);
                // b*c
                single_multiply(hstackpos+2*i+1, lstackpos+2*i+1, hstackpos[2*i+old_len], lstackpos[2*i+old_len]);
                // b*d
                single_multiply(&bh, &bl, hstackpos[2*i+old_len+1], lstackpos[2*i+old_len+1]);
                // a*d
                single_multiply(&ah, &al, hstackpos[2*i+old_len+1], lstackpos[2*i+old_len+1]);
                hstackpos[2*i] -= bl;
                lstackpos[2*i] -= bh;
                hstackpos[2*i+1] += ah;
                lstackpos[2*i+1] += al;
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
        int k=old_len+len*step;
        for (int i=len-1; i >= 0; i--) {
            for (int j=step-1; j >= 0; j--) {
                k--;
                hstackpos[k] = hstackpos[old_len+i];
                lstackpos[k] = lstackpos[old_len+i];
                single_multiply(hstackpos+k, lstackpos+k, hstackpos[j], lstackpos[j]);   //pos[i*step+j] = pos[i]*stackpos[j];
            }
        }
        *result_length = len;
        *result_type = type | *result_type;
        for (int i=0; i < len*step; i++) hstackpos[i] = hstackpos[i+old_len];
        for (int i=0; i < len*step; i++) lstackpos[i] = lstackpos[i+old_len];
    } else {
        // Otherwise, we broadcast and add. We also know that if type is a
        // list, then result_type is as well, otherwise we would be in the
        // previous if-block.
        if (IS_TYPE(type, TYPE_POINT)) {
            // If type is a point, then result_type is not a point
            if (type & TYPE_LIST) {
                // List of scalars times a list of points. result_length has already
                // been shortened to the appropriate length. The case where the result
                // is not a list has already been handled above.
                for (int i=0; i < *result_length; i++) {
                    single_multiply(hstackpos+old_len+2*i, lstackpos+old_len+2*i, hstackpos[i], lstackpos[i]);
                    single_multiply(hstackpos+old_len+2*i+1, lstackpos+old_len+2*i+1, hstackpos[i], lstackpos[i]);
                }
                memcpy(hstackpos, hstackpos+old_len, (*result_length)*step*sizeof(double));
                memcpy(lstackpos, lstackpos+old_len, (*result_length)*step*sizeof(double));
            } else {
                // Scalar or list of scalars times a single point
                double xh = hstackpos[old_len], xl = lstackpos[old_len], yh = hstackpos[old_len+1], yl = lstackpos[old_len+1];
                for (int i=*result_length-1; i>=0; i--) {
                    // y comes first because x might overwrite stackpos[0]
                    hstackpos[2*i+1] = hstackpos[i]; hstackpos[2*i] = hstackpos[i];
                    lstackpos[2*i+1] = lstackpos[i]; lstackpos[2*i] = lstackpos[i];
                    single_multiply(hstackpos+2*i, lstackpos+2*i, xh, xl);
                    single_multiply(hstackpos+2*i+1, lstackpos+2*i+1, yh, yl);
                }
            }
            *result_type |= TYPE_POINT;
        } else {
            for (int i=0; i < *result_length; i++) {
                for (int j=0; j < step; j++) single_multiply(hstackpos+i*step+j, lstackpos+i*step+j, hstackpos[old_len+(i%len)], lstackpos[old_len+(i%len)]);
            }
        }
    }
}

uint32_t interval_multiply(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t type;
    uint32_t len;
    uint32_t result_type = 0;
    uint32_t result_length = 0;
    int i=0;
    result_type = arg->inter(arg, hstackpos, lstackpos);
    result_length = result_type>>8;
    result_length /= GET_STEP(result_type);
    arg = arg->next_arg;
    //printf("Initial value (%08x): ", result_type); print_object(result_type, stackpos); printf("\n");
    while (arg) {
        type = arg->inter(arg, hstackpos+result_length*GET_STEP(result_type), lstackpos+result_length*GET_STEP(result_type));
        // Keep track of the position in case result_length is changed.
        interval_multiply_in_place(hstackpos, lstackpos, &result_type, &result_length, type);
        arg = arg->next_arg;
    }
    result_length *= GET_STEP(result_type);
    apply_sign(fs->value_type, hstackpos, lstackpos, 0, result_length);
    
    return (result_length<<8) | (result_type&0xff);
}

void mipow(double *hstackpos, double *lstackpos, double h, double l) {
    // pow(a, b)
    // This is defined when a>0 or when a=0 and b>=0 or when a<0 and b is an integer
    // Intervals must only contain the values for which the function is defined.
    // pow(a, b) can be negative when b is odd.
    double al = lstackpos[0], ah = hstackpos[0];
    double ih = floor(h), il = ceil(l);
    double ahbh = pow(ah, h), ahbl = pow(ah, l), albh = pow(al, h), albl = pow(al, l);
    if (al < 0) {
        albh = pow(al, ih);
        albl = pow(al, il);
    }
    if (ah < 0) {
        ahbh = pow(ah, ih);
        ahbl = pow(ah, il);
    }

    lstackpos[0] = ahbh;
    hstackpos[0] = ahbh;
    if ((ah >= 0) && (al <= 0)) lstackpos[0] = 0;
    if (ahbl > hstackpos[0]) hstackpos[0] = ahbl;
    if (ahbl < lstackpos[0]) lstackpos[0] = ahbl;
    if (albh > hstackpos[0]) hstackpos[0] = albh;
    if (albh < lstackpos[0]) lstackpos[0] = albh;
    if (albl > hstackpos[0]) hstackpos[0] = albl;
    if (albl < lstackpos[0]) lstackpos[0] = albl;
}

uint32_t interval_exponentiate(void *f, double *hstackpos, double *lstackpos) {
    return interval_general_two_args(f, hstackpos, lstackpos, mipow);
}

uint32_t interval_div(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t t1, t2, l1, l2;
    uint8_t step1, step2;
    t1 = arg->inter(arg, hstackpos, lstackpos);
    l1 = t1>>8;
    step1 = GET_STEP(t1);

    arg = arg->next_arg;
    t2 = arg->inter(arg, hstackpos+l1, lstackpos+l1);
    l2 = t2>>8;
    step2 = GET_STEP(t2);

    uint32_t result_len = l1/step1;
    double val;
    if (((l2/step2 < result_len) || !(t1 & TYPE_LIST)) && (t2 & TYPE_LIST)) result_len = l2/step2;
    if (IS_TYPE(t1, TYPE_DOUBLE)) {
        if (IS_TYPE(t2, TYPE_DOUBLE)) {
            // Dividing real numbers
            if (!(t1 & TYPE_LIST)) {
                double vh = hstackpos[0];
                double vl = lstackpos[0];
                for (int i=0; i < result_len; i++) {
                    if ((hstackpos[l1+i] >= 0) && (lstackpos[l1+i] <= 0)) {
                        hstackpos[i] = INFINITY;
                        lstackpos[i] = -INFINITY;
                    } else {
                        hstackpos[i] = vh;
                        lstackpos[i] = vl;
                        single_multiply(hstackpos+i, lstackpos+i, 1/lstackpos[l1+i], 1/hstackpos[l1+i]);
                    }
                }
                apply_sign(fs->value_type, hstackpos, lstackpos, 0, result_len);
                return (result_len<<8) | ((t1 | t2) & 0xff);
            }
            for (int i=0; i < result_len; i++) {
                if ((hstackpos[l1+i] >= 0) && (lstackpos[l1+i] <= 0)) {
                    hstackpos[i] = INFINITY;
                    lstackpos[i] = -INFINITY;
                } else single_multiply(hstackpos+i, lstackpos+i, 1/lstackpos[l1+i], 1/hstackpos[l1+i]);
            }
            apply_sign(fs->value_type, hstackpos, lstackpos, 0, result_len);
            return (result_len<<8) | ((t1 | t2) & 0xff);
        } /*else {
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
        return (result_len<<9) | ((t1 | t2) & 0xff);*/
    }
    return 0;
}

void miequals(double *hstackpos, double *lstackpos, double h, double l) {
    hstackpos[0] -= l;
    lstackpos[0] -= h;
}

uint32_t interval_equals(void *f, double *hstackpos, double *lstackpos) {
    return interval_general_two_args(f, hstackpos, lstackpos, miequals);
}

void migt(double *hstackpos, double *lstackpos, double h, double l) {
    hstackpos[0] -= l;
    lstackpos[0] -= h;
}

void milt(double *hstackpos, double *lstackpos, double h, double l) {
    l -= hstackpos[0];
    h -= lstackpos[0];
    hstackpos[0] = h;
    lstackpos[0] = l;
}

uint32_t interval_compare_single(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint64_t cmp = *((uint64_t*)(fs->value));
    if (cmp & 1) return interval_general_two_args(f, hstackpos, lstackpos, milt);
    else return interval_general_two_args(f, hstackpos, lstackpos, migt);
}

uint32_t interval_compare(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint64_t cmp = *((uint64_t*)(fs->value));
    uint32_t result_type = 0;
    uint32_t result_len = 0;
    uint32_t last_len = 0;
    uint32_t type, len, newlen;
    result_type = arg->inter(arg, hstackpos, lstackpos);
    len = result_type>>8;
    memcpy(hstackpos+len, hstackpos, len*sizeof(double));
    memcpy(lstackpos+len, lstackpos, len*sizeof(double));
    last_len = len; result_len = len;
    arg = arg->next_arg;
    while (arg) {
        type = arg->inter(arg, hstackpos+result_len+last_len, lstackpos+result_len+last_len);
        len = type>>8;
        if ((type & TYPE_LIST) && !(result_type & TYPE_LIST)) {
            // If the result is not a list but the operand is, expand the result array
            for (int i=last_len+len-1; i >= 0; i--) hstackpos[len+i] = hstackpos[result_len+i];
            for (int i=last_len+len-1; i >= 0; i--) lstackpos[len+i] = lstackpos[result_len+i];
            for (int i=1; i < len; i++) hstackpos[i] = hstackpos[0];
            for (int i=1; i < len; i++) lstackpos[i] = lstackpos[0];
            result_len = len;
            result_type |= TYPE_LIST;
        }
        newlen = result_len;
        if (!(result_type & TYPE_LIST) || ((type & TYPE_LIST) && (len < result_len))) newlen = len;
        int i2, i3;
        double temp;
        if (cmp & 1) {
            // less than
            for (int i=newlen-1; i >= 0; i--) {
                i2 = result_len+i%last_len;
                i3 = result_len+last_len+i%len;
                temp = hstackpos[i2];
                hstackpos[i2] = hstackpos[i3] - lstackpos[i2];
                lstackpos[i2] = lstackpos[i3] - temp;
                if (hstackpos[i2] < hstackpos[i]) hstackpos[i] = hstackpos[i2];
                if (lstackpos[i2] < lstackpos[i]) lstackpos[i] = lstackpos[i2];
            }
        } else {
            // greater than
            for (int i=newlen-1; i >= 0; i--) {
                i2 = result_len+i%last_len;
                i3 = result_len+last_len+i%len;
                hstackpos[i2] -= lstackpos[i3];
                lstackpos[i2] -= hstackpos[i3];
                if (hstackpos[i2] < hstackpos[i]) hstackpos[i] = hstackpos[i2];
                if (lstackpos[i2] < lstackpos[i]) lstackpos[i] = lstackpos[i2];
            }
        }
        cmp = cmp >> 2;
        for (int i=0; i < len; i++) hstackpos[newlen+i] = hstackpos[result_len+last_len+i];
        for (int i=0; i < len; i++) lstackpos[newlen+i] = lstackpos[result_len+last_len+i];
        result_len = newlen;
        last_len = len;
        arg = arg->next_arg;
    }
    return (result_len << 8) | (result_type & TYPE_LIST) | TYPE_BOOLEAN;
}


uint32_t interval_greater(void *f, double *hstackpos, double *lstackpos) {
    return interval_general_two_args(f, hstackpos, lstackpos, migt);
}

void micos(double *hstackpos, double *lstackpos) {
    double h = hstackpos[0], l = lstackpos[0];
    if (h - l >= 2*M_PI) {
        hstackpos[0] = 1;
        lstackpos[0] = -1;
    } else {
        double hcos = cos(h), lcos = cos(l);
        if (hcos > lcos) {
            hstackpos[0] = hcos;
            lstackpos[0] = lcos;
        } else {
            hstackpos[0] = lcos;
            lstackpos[0] = hcos;
        }
        // cos(x*M_PI) is 1 when x is even and -1 when x is odd
        h = h / M_PI;
        int32_t l_int = ceil(l / M_PI);
        if (l_int % 2) {
            // an odd number must be in range, so set the minimum
            if (l_int < h) lstackpos[0] = -1;
            // If an even number is also in range, set the maximum
            if (l_int+1 < h) hstackpos[0] = 1;
        } else {
            if (l_int < h) hstackpos[0] = 1;
            if (l_int+1 < h) lstackpos[0] = -1;
        }
    }
}

void misin(double *hstackpos, double *lstackpos) {
    hstackpos[0] -= M_PI_2;
    lstackpos[0] -= M_PI_2;
    micos(hstackpos, lstackpos);
}

void mitan(double *hstackpos, double *lstackpos) {
    if (ceil(lstackpos[0]/M_PI - 0.5) <= floor(hstackpos[0]/M_PI - 0.5)) {
        hstackpos[0] = INFINITY;
        lstackpos[0] = -INFINITY;
    } else {
        hstackpos[0] = tan(hstackpos[0]);
        lstackpos[0] = tan(lstackpos[0]);
    }
}

void miarctan(double *hstackpos, double *lstackpos) {
    hstackpos[0] = atan(hstackpos[0]);
    lstackpos[0] = atan(lstackpos[0]);
}

void miexp(double *hstackpos, double *lstackpos) {
    hstackpos[0] = exp(hstackpos[0]);
    lstackpos[0] = exp(lstackpos[0]);
}

uint32_t interval_cosine(void *f, double *hstackpos, double *lstackpos) {
    return interval_general_one_arg(f, hstackpos, lstackpos, micos);
}

uint32_t interval_sine(void *f, double *hstackpos, double *lstackpos) {
    return interval_general_one_arg(f, hstackpos, lstackpos, misin);
}

uint32_t interval_tangent(void *f, double *hstackpos, double *lstackpos) {
    return interval_general_one_arg(f, hstackpos, lstackpos, mitan);
}

uint32_t interval_arctan(void *f, double *hstackpos, double *lstackpos) {
    return interval_general_one_arg(f, hstackpos, lstackpos, miarctan);
}

uint32_t interval_exp(void *f, double *hstackpos, double *lstackpos) {
    return interval_general_one_arg(f, hstackpos, lstackpos, miexp);
}

uint32_t interval_list(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    if (!arg) return TYPE_LIST;
    uint8_t type = 0xff;
    uint32_t argtype;
    uint32_t arglen;
    uint32_t st=0;
    uint8_t ellipsis = 0;
    /*do {
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
    } while (arg);*/
    while (arg) {
        argtype = arg->inter(arg, hstackpos+st, lstackpos+st);
        if (argtype & TYPE_LIST) FAIL("ERROR: Cannot store a list inside a list\n");
        if ((type != 0xff) && ((argtype&0xff) != type)) FAIL("ERROR: All list elements must have the same type\n");
        st += argtype>>8;
        type = argtype & 0xff;
        arg = arg->next_arg;
    }
    return (st<<8) | type | TYPE_LIST;
}

uint32_t interval_point(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    uint32_t t1, t2, l1, l2;
    t1 = arg->inter(arg, hstackpos, lstackpos);
    l1 = t1 >> 8;
    arg = arg->next_arg;
    t2 = arg->inter(arg, hstackpos+l1, lstackpos+l1);
    l2 = t2>>8;
    
    // Faster for parametric functions, no memory allocation needed
    if (!(t1 & TYPE_LIST) && !(t2 & TYPE_LIST)) {
        //printf("Interval of point, ([%f, %f], [%f, %f])\n", lstackpos[0], hstackpos[0], lstackpos[1], hstackpos[1]);
        apply_sign(fs->value_type, hstackpos, lstackpos, 0, 2);
        return (2<<8) | TYPE_POINT;
    }

    double *htemp = malloc(2*(l1+l2)*sizeof(double));
    double *ltemp = htemp + l1+l2;
    memcpy(htemp, hstackpos, (l1+l2)*sizeof(double));
    memcpy(ltemp, lstackpos, (l1+l2)*sizeof(double));
    
    uint32_t result_len = l1;
    if (!(t1 & TYPE_LIST) || ((t2 & TYPE_LIST) && (l2 < result_len))) result_len = l2;
    for (int i=0; i < result_len; i++) {
        hstackpos[2*i] = htemp[i%l1];
        hstackpos[2*i+1] = htemp[l1+i%l2];
    }
    for (int i=0; i < result_len; i++) {
        lstackpos[2*i] = ltemp[i%l1];
        lstackpos[2*i+1] = ltemp[l1+i%l2];
    }
    apply_sign(fs->value_type, hstackpos, lstackpos, 0, 2*result_len);
    return (result_len<<9) | TYPE_POINT;
}

uint32_t interval_extract_x(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    //printf("extracting x, arg %p, oper %p (%p), value %p, value_type %08x\n", arg, arg->oper, (func_value), arg->value, arg->value_type);
    
    uint32_t argtype, arglen;
    argtype = arg->inter(arg, hstackpos, lstackpos);
    arglen = argtype>>8;
    //printf("extracting from type %08x\n", argtype);
    if (!((argtype & TYPE_MASK) == TYPE_POINT)) FAIL("ERROR: cannot access x-coordinate of type %08x, block %p\n", argtype, fs);
    for (int i=0; i < arglen; i += 2) hstackpos[i>>1] = hstackpos[i];
    for (int i=0; i < arglen; i += 2) lstackpos[i>>1] = lstackpos[i];
    return ((argtype>>1) & 0xffffff00) | ((argtype&0xff) & ~TYPE_MASK);
}

uint32_t interval_extract_y(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    
    uint32_t argtype, arglen;
    argtype = arg->inter(arg, hstackpos, lstackpos);
    arglen = argtype>>8;
    if (!((argtype & TYPE_MASK) == TYPE_POINT)) FAIL("ERROR: cannot access y-coordinate of type %08x\n", argtype);
    for (int i=0; i < arglen; i += 2) hstackpos[i>>1] = hstackpos[i+1];
    for (int i=0; i < arglen; i += 2) lstackpos[i>>1] = lstackpos[i+1];
    return ((argtype>>1) & 0xffffff00) | ((argtype&0xff) & ~TYPE_MASK);
}

uint32_t interval_user_defined(void *f, double *hstackpos, double *lstackpos) {
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
        type = arg->inter(arg, hstackpos, lstackpos+st);
        len = type>>8;
        memcpy(lstackpos+st+len, hstackpos, len*sizeof(double));
        //printf("argument is %p [%f, %f]\n", arg, lstackpos[st], hstackpos[st]);
        target_arg[varnum].pointer = lstackpos+st;
        target_arg[varnum].type = type;
        target_arg[varnum].flags |= VARIABLE_INTERVAL;
        st += 2*len;
        //printf("argument is %08x:", type); print_object(type, target_arg[varnum].pointer); printf(" to "); print_object(type, target_arg[varnum].pointer+len); printf("\n");
        varnum++;
        arg = arg->next_arg;
    } while (arg);
    type = target->inter(target, hstackpos+st, lstackpos+st);
    //printf("target %p, type %08x, [%f, %f]\n", target, type, lstackpos[st], hstackpos[st]);
    len = type>>8;
    // Shrink the stack by removing unused values
    apply_sign(fs->value_type, hstackpos, lstackpos, st, len);
    //printf("return type is %08x\n", type); print_object(type, stackpos); printf("\n");
    return type;
}

double integrate_func_interval(double x, void *params) {
    function *expr = (function*)(((void**)params)[0]);
    double *hstackpos = (double*)(((void**)params)[1]);
    double *lstackpos = (double*)(((void**)params)[2]);
    uint8_t *flags = (uint8_t*)(((void**)params)[4]);
    ((variable*)(((void**)params)[3]))->pointer = &x;
    expr->inter(expr, hstackpos, lstackpos);
    //printf("interval of integrant is [%f, %f] at %f\n", lstackpos[0], hstackpos[0], x);
    if (hstackpos[0] != lstackpos[0]) *flags &= 0x7f;
    if (*flags & 0x01) return hstackpos[0];
    else return lstackpos[0];
}

uint32_t interval_integrate_gsl(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *var = fs->first_arg;
    function *lb = var->next_arg;
    function *ub = lb->next_arg;
    function *expr = ub->next_arg;

    lb->inter(lb, hstackpos, lstackpos);
    double lbv[2] = {lstackpos[0], hstackpos[0]};
    ub->inter(ub, hstackpos, lstackpos);
    double ubv[2] = {lstackpos[0], hstackpos[0]};
    variable *varp = (variable*)(var->value);
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    uint8_t flags = 0x80;
    ((variable*)(var->value))->type = 1<<8;
    void *params[5] = {expr, hstackpos, lstackpos, var->value, &flags};
    F.function = &integrate_func_interval;
    F.params = params;
    double hresult, lresult, error;
    int hstatus=0, lstatus=0;
    gsl_set_error_handler_off();
    double diffh, diffl;
    // No overlap, evaluate on the smallest possible range. If this is not valid,
    // then the interval will be [NaN, NaN], indicating that the function is defined
    // nowhere on the given intervals
    double min_lbv = lbv[1];
    double min_ubv = ubv[0];
    double lstar, ustar;
    double temp[2];
    if (ubv[0] >= lbv[1]) {
        lstar = lbv[1];
        ustar = ubv[0];
    } else if (lbv[0] >= ubv[1]) {
        lstar = lbv[0];
        ustar = ubv[1];
    } else {
        // lbv and ubv intervals overlap
        lstar = (fmin(ubv[1], lbv[1]) + fmax(ubv[0], ubv[0]))/2;
        ustar = lstar;
    }
    lstatus = gsl_integration_qags(&F, lstar, ustar, 0, 1e-7, 1000, w, &lresult, &error);
    if (!(flags & 0x80)) {
        //printf("computing hresult separately\n");
        flags |= 0x01;
        hstatus = gsl_integration_qags(&F, lstar, ustar, 0, 1e-7, 1000, w, &hresult, &error);
    } else hresult = lresult;
    //printf("lresult %f, hresult %f, lstar %f, ustar %f\n", lresult, hresult, lstar, ustar);
    if (lresult > hresult) {
        temp[0] = hresult;
        hresult = lresult;
        lresult = temp[0];
    }
    if (lstatus || hstatus) {
        // Minimum interval is invalid
        hstackpos[0] = NAN;
        lstackpos[0] = NAN;
        return 1<<8;
    }

    varp->pointer = temp;
    varp->flags |= VARIABLE_INTERVAL;
    diffh = 0; diffl = 0;
    if (lbv[0] != lstar) {
        temp[0] = lbv[0]; temp[1] = lstar;
        expr->inter(expr, hstackpos, lstackpos);
        diffh += lstackpos[0]*(lbv[0] - lstar);
        diffl += hstackpos[0]*(lbv[0] - lstar);
    }
    if (lstar != lbv[1]) {
        temp[0] = lstar; temp[1] = lbv[1];
        expr->inter(expr, hstackpos, lstackpos);
        lstackpos[0] *= (lbv[1] - lstar);
        hstackpos[0] *= (lbv[1] - lstar);
        if (lstackpos[0] < diffl) diffl = lstackpos[0];
        if (hstackpos[0] > diffh) diffh = hstackpos[0];
    }
    lresult += diffl;
    hresult += diffh;
    diffl = 0; diffh = 0;
    if (ubv[0] != ustar) {
        temp[0] = ubv[0]; temp[1] = ustar;
        expr->inter(expr, hstackpos, lstackpos);
        diffh += lstackpos[0]*(ubv[0] - ustar);
        diffl += hstackpos[0]*(ubv[0] - ustar);
    }
    if (ustar != ubv[1]) {
        temp[0] = ustar; temp[1] = ubv[1];
        expr->inter(expr, hstackpos, lstackpos);
        lstackpos[0] *= (ubv[1] - ustar);
        hstackpos[0] *= (ubv[1] - ustar);
        if (lstackpos[0] < diffl) diffl = lstackpos[0];
        if (hstackpos[0] > diffh) diffh = hstackpos[0];
    }
    lresult += diffl;
    hresult += diffh;
    varp->flags &= ~VARIABLE_INTERVAL;
    hstackpos[0] = hresult;
    lstackpos[0] = lresult;
    return 1<<8;
}

