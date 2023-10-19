#include "linalg_functions.h"
#include "functions.h"
#include "parse.h"
#include <stdint.h>
#include <stdio.h>
#include <lapack.h>

#define SIGN_BIT(v) (((v->value_type)&0x80) ? -1 : 1)

#define FAIL(...) {printf(__VA_ARGS__); exit(EXIT_FAILURE);}
#define IS_TYPE(type, ref) (((type) & TYPE_MASK) == (ref))

uint32_t func_solve(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    printf("solving system of equations, arg %p\n", arg);
    
    uint32_t st=0;
    uint32_t argtype, alen, blen;
    int nrhs = 1;
    double *a = stackpos;
    argtype = arg->oper(arg, stackpos);
    alen = argtype>>8;
    st += alen;
    if (!IS_TYPE(argtype, TYPE_DOUBLE)) FAIL("ERROR: A matrix must be array of doubles\n");
    printf("first argument: "); print_object(argtype, stackpos);

    arg = arg->next_arg;
    double *b = stackpos+st;
    argtype = arg->oper(arg, stackpos+st);
    blen = argtype>>8;
    st += blen;
    if (!IS_TYPE(argtype, TYPE_DOUBLE)) FAIL("ERROR: B matrix must be array of doubles\n");
    if (arg->next_arg) {
        arg = arg->next_arg;
        argtype = arg->oper(arg, stackpos+st);
        nrhs = (int)stackpos[st];
    }
    if (blen % nrhs) FAIL("ERROR: Bad dimension for B matrix, nrhs %d, blen %d\n", nrhs, blen);

    int arows = blen / nrhs;
    int acols = alen / arows;
    int ldb = arows;
    printf("nrhs %d, arows %d, acols %d, ldb %d\n", nrhs, arows, acols, ldb);
    if (acols > ldb) {
        ldb = acols;
        // Space out the columns of B to make room for the result
        for (int i=nrhs-1; i>=0; i--) {
            for (int j=arows-1; j>=0; j--) {
                b[i*ldb + j] = b[i*arows + j];
            }
        }
        st = alen + ldb*nrhs;
    }
    int *jpvt = (int*)(stackpos+st);
    st += acols;
    
    int lwork = -1;
    double rcond = 1e-14;
    int rank, info;
    LAPACK_dgelsy(&arows, &acols, &nrhs, a, &arows, b, &ldb, jpvt, &rcond, &rank, stackpos+st, &lwork, &info);
    lwork = stackpos[st];
    LAPACK_dgelsy(&arows, &acols, &nrhs, a, &arows, b, &ldb, jpvt, &rcond, &rank, stackpos+st, &lwork, &info);
    for (int i=0; i < nrhs; i++) {
        for (int j=0; j < acols; j++) {
            stackpos[i*acols + j] = b[i*ldb + j]*SIGN_BIT(fs);
        }
    }
    return (nrhs*acols)<<8;
}
