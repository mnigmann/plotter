#include "linalg_intervals.h"
#include "linalg_functions.h"
#include "functions.h"
#include "intervals.h"
#include "parse.h"
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <lapack.h>
#include <complex.h>
#include <string.h>

#define FAIL(...) {printf(__VA_ARGS__); exit(EXIT_FAILURE);}
#define IS_TYPE(type, ref) (((type) & TYPE_MASK) == (ref))

uint32_t interval_det(void *f, double *hstackpos, double *lstackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;

    uint32_t argtype, alen, blen;
    argtype = arg->inter(arg, hstackpos, lstackpos);
    alen = argtype>>8;
    int info;
    if (!IS_TYPE(argtype, TYPE_DOUBLE)) FAIL("ERROR: interval_det only supports real matrices\n");

    int32_t n = isqrt(alen);
    double *work = hstackpos+alen;
    int *ipiv = (int*)(work+alen);

    if (n*n != alen) FAIL("ERROR: can only compute determinant of square matrix\n");

    // Compute sum and difference matrices
    double temp;
    for (int i=0; i < alen; i++) {
        temp = hstackpos[i];
        hstackpos[i] = (hstackpos[i]+lstackpos[i])/2;
        lstackpos[i] = (temp - lstackpos[i])/2;
    }

    LAPACK_dgetrf(&n, &n, hstackpos, &n, ipiv, &info);
    double det = 1;
    for (int i=0; i < n; i++) {
        det *= hstackpos[(n+1)*i];
        if (i+1 != ipiv[i]) det = -det;
    }
    LAPACK_dgetri(&n, hstackpos, &n, ipiv, work, (int32_t*)(&alen), &info);
    
    // Decompose outer product
    // Compute maximum in each row. This is the u vector
    for (int i=0; i < n; i++) {
        work[i] = lstackpos[i*n];
        work[i+n] = 0;
        for (int j=1; j < n; j++) {
            if (lstackpos[i*n+j] > work[i]) work[i] = lstackpos[i*n+j];
        }
    }
    // Compute v vector
    for (int i=0; i < n; i++) {
        for (int j=0; j < n; j++) {
            if (lstackpos[i*n+j]/work[i] > work[j+n]) work[j+n] = lstackpos[i*n+j]/work[i];
        }
    }
    // Evaluate v @ |inv((H+L)/2)| @ u
    double prod = 0;
    for (int i=0; i < n; i++) {
        for (int j=0; j < n; j++) {
            prod += fabs(hstackpos[i*n+j])*work[j]*work[i+n];
        }
    }
    // Apply matrix determinant lemma
    double d1 = det*(1+prod);
    double d2 = det*(1-prod);
    if (d1 > d2) {
        hstackpos[0] = d1;
        lstackpos[0] = d2;
    } else {
        hstackpos[0] = d2;
        lstackpos[0] = d1;
    }
    //printf("estimate is %f, interval [%f, %f]\n", prod, lstackpos[0], hstackpos[0]);
    return 1<<8;
}

/*int main() {
    double hmat[] = {10,  2,  2,  9,  3,  
                      7, 10, 10,  8,  7,  
                      9,  2,  3, 11,  5,  
                      4,  2,  7,  4,  7, 
                     10,  7,  5,  2,  7};
    double lmat[] = { 6,  0,  0,  9,  1,  
                      7,  6,  8,  6,  5,  
                      9,  2, -1,  7,  3,  
                      0,  2,  7,  0,  5,  
                      6,  5,  1, -2,  7};
    
    double hstack[1024];
    double lstack[1024];
    memcpy(hstack, lmat, 25*sizeof(double));
    memcpy(hstack+25, hmat, 25*sizeof(double));

    variable v = new_variable("var", (25<<8), VARIABLE_INTERVAL, hstack);
    function vb = new_function(NULL, NULL, NULL);
    vb.inter = interval_value;
    vb.value_type = 0x40;
    vb.value = (void*)(&v);
    function detb = new_function(NULL, NULL, &vb);
    detb.inter = interval_det;

    detb.inter(&detb, hstack+25, lstack+25);
}*/
