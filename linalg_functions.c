#include "linalg_functions.h"
#include "functions.h"
#include "parse.h"
#include <stdint.h>
#include <stdio.h>
#include <lapack.h>
#include <complex.h>
#include <string.h>

#define SIGN_BIT(v) (((v->value_type)&0x80) ? -1 : 1)

#define FAIL(...) {printf(__VA_ARGS__); exit(EXIT_FAILURE);}
#define IS_TYPE(type, ref) (((type) & TYPE_MASK) == (ref))
#define CSPLIT(z) creal(z),cimag(z)

uint32_t func_solve(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    //printf("solving system of equations, arg %p\n", arg);
    
    uint32_t st=0;
    uint32_t argtype, alen, blen;
    int nrhs = 1;
    double *a = stackpos;
    argtype = arg->oper(arg, stackpos);
    alen = argtype>>8;
    st += alen;
    if (!IS_TYPE(argtype, TYPE_DOUBLE)) FAIL("ERROR: A matrix must be array of doubles\n");
    //printf("first argument: "); print_object(argtype, stackpos); printf("\n");

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

    if ((alen == 0) || (blen == 0)) return TYPE_LIST;

    int arows = blen / nrhs;
    int acols = alen / arows;
    int ldb = arows;
    //printf("nrhs %d, arows %d, acols %d, ldb %d\n", nrhs, arows, acols, ldb);
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
    return ((nrhs*acols)<<8) | TYPE_LIST;
}

uint32_t isqrt(uint32_t n) {
    uint32_t x = n;
    uint32_t c = 0;
    uint32_t d = 1<<30;
    while (d > n) d>>=2;
    while (d != 0) {
        if (x >= c+d) {
            x -= c+d;
            c = (c>>1) + d;
        } else {
            c >>= 1;
        }
        d >>= 2;
    }
    return c;
}

void compute_eigvals(double *stackpos, int32_t n) {
    uint32_t alen = 2*n*n;
    uint32_t blen = alen;
    lapack_int info, lwork=-1, lda=n, ldb=n, ldv=1;
    lapack_complex_double *result = (lapack_complex_double*)(stackpos+alen+blen);
    lapack_complex_double lwork_out = -1;
    lapack_complex_double *alpha = malloc(n*sizeof(lapack_complex_double));
    lapack_complex_double *beta = malloc(n*sizeof(lapack_complex_double));
    LAPACK_zggev("N", "N", &n, (lapack_complex_double*)stackpos, &n, (lapack_complex_double*)(stackpos+alen), &n, result, result+n, NULL, &ldv, NULL, &ldv, &lwork_out, &lwork, (double*)(result+2*n), &info);
    lwork = creal(lwork_out);
    lapack_complex_double *work = malloc(lwork*sizeof(lapack_complex_double));
    LAPACK_zggev("N", "N", &n, (lapack_complex_double*)stackpos, &n, (lapack_complex_double*)(stackpos+alen), &n, alpha, beta, NULL, &ldv, NULL, &ldv, work, &lwork, (double*)(result+2*n), &info);
    free(work);
    if (info != 0) printf("WARNING: zggev produced non-zero error code %d\n", info);
    lapack_complex_double lambda;
    for (int i=0; i < n; i++) {
        //lambda = result[i]/result[i+n];
        lambda = alpha[i]/beta[i];
        //printf("eigval is %e%+ei, %e%+ei, %e%+ei\n", CSPLIT(alpha[i]), CSPLIT(beta[i]), CSPLIT(lambda));
        stackpos[2*i] = creal(lambda);
        stackpos[2*i+1] = cimag(lambda);
    }
    free(alpha);
    free(beta);
}

void compute_eigvals_real(double *stackpos, int32_t n) {
    uint32_t alen = 2*n*n;
    uint32_t blen = alen;
    for (int i=0; i < alen; i+=2) stackpos[i/2] = stackpos[i];
    for (int i=0; i < blen; i+=2) stackpos[alen+i/2] = stackpos[alen+i];
    lapack_int info, lwork=-1, lda=n, ldb=n, ldv=1;
    double lwork_out = -1;
    double *alpha_r = malloc(n*sizeof(double));
    double *alpha_i = malloc(n*sizeof(double));
    double *beta = malloc(n*sizeof(double));
    LAPACK_dggev("N", "N", &n, stackpos, &n, stackpos+alen, &n, alpha_r, alpha_i, beta, NULL, &ldv, NULL, &ldv, &lwork_out, &lwork, &info);
    lwork = creal(lwork_out);
    printf("lwork is %d, info %d\n", lwork, info);
    double *work = malloc(lwork*sizeof(double));
    LAPACK_dggev("N", "N", &n, stackpos, &n, stackpos+alen, &n, alpha_r, alpha_i, beta, NULL, &ldv, NULL, &ldv, work, &lwork, &info);
    free(work);
    if (info != 0) printf("WARNING: zggev produced non-zero error code %d\n", info);
    lapack_complex_double lambda;
    for (int i=0; i < n; i++) {
        //lambda = result[i]/result[i+n];
        lambda = (alpha_r[i] + I*alpha_i[i])/beta[i];
        //printf("eigval is %e%+ei, %e, %e%+ei\n", alpha_r[i], alpha_i[i], beta[i], CSPLIT(lambda));
        stackpos[2*i] = creal(lambda);
        stackpos[2*i+1] = cimag(lambda);
    }
    free(alpha_i);
    free(alpha_r);
    free(beta);
}

uint32_t func_eigvals(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;
    
    uint32_t argtype, alen, blen;
    argtype = arg->oper(arg, stackpos);
    alen = argtype>>8;
    if (!IS_TYPE(argtype, TYPE_POINT)) {
        // Real matrix provided. Convert to complex
        for (int i=alen-1; i>=0; i--) {
            stackpos[2*i+1] = 0;
            stackpos[2*i] = stackpos[i];
        }
        alen = alen * 2;
    }
    int32_t n = isqrt(alen / 2);
    if (n*n*2 != alen) FAIL("ERROR: matrix must be square to compute eigenvalues\n");
    arg = arg->next_arg;
    int info, lwork=-1, lda=n, ldb=n, ldv=1;
    if (arg) {
        // Two matrices provided, compute generalized eigenvalues
        argtype = arg->oper(arg, stackpos+alen);
        blen = argtype>>8;
        if (!IS_TYPE(argtype, TYPE_POINT)) {
            // Real matrix provided. Convert to complex
            for (int i=blen-1; i>=0; i--) {
                stackpos[alen+2*i+1] = 0;
                stackpos[alen+2*i] = stackpos[alen+i];
            }
            blen = blen * 2;
        }
        if (alen != blen) FAIL("ERROR: matrices must have the same dimension to compute generalized eigenvalues\n");
        memset(stackpos+alen+blen, 0, sizeof(double)*12*n);
        
        compute_eigvals_real(stackpos, n);

        /*lapack_complex_double *result = (lapack_complex_double*)(stackpos+alen+blen);
        lapack_complex_double lwork_out = -1;
        lapack_complex_double *alpha = malloc(n*sizeof(lapack_complex_double));
        lapack_complex_double *beta = malloc(n*sizeof(lapack_complex_double));
        LAPACK_zggev("N", "N", &n, (lapack_complex_double*)stackpos, &n, (lapack_complex_double*)(stackpos+alen), &n, result, result+n, NULL, &ldv, NULL, &ldv, &lwork_out, &lwork, (double*)(result+2*n), &info);
        lwork = creal(lwork_out);
        printf("lwork is %d, info %d\n", lwork, info);
        lapack_complex_double *work = malloc(lwork*sizeof(lapack_complex_double));
        LAPACK_zggev("N", "N", &n, (lapack_complex_double*)stackpos, &n, (lapack_complex_double*)(stackpos+alen), &n, alpha, beta, NULL, &ldv, NULL, &ldv, work, &lwork, (double*)(result+2*n), &info);
        free(work);
        if (info != 0) printf("WARNING: zggev produced non-zero error code %d\n", info);
        for (int i=0; i < alen+blen+12*n+lwork*2; i++) {
            printf("%e\t", stackpos[i]);
            if (i%8 == 7) printf("\n");
        }
        lapack_complex_double lambda;
        for (int i=0; i < n; i++) {
            //lambda = result[i]/result[i+n];
            lambda = alpha[i]/beta[i];
            printf("eigval is %e%+ei, %e%+ei, %e%+ei\n", CSPLIT(alpha[i]), CSPLIT(beta[i]), CSPLIT(lambda));
            stackpos[2*i] = creal(lambda);
            stackpos[2*i+1] = cimag(lambda);
        }
        free(alpha);
        free(beta);*/
        for (int i=0; i < n*2; i++) stackpos[i] *= SIGN_BIT(fs);

        return (n<<9) | TYPE_POINT | TYPE_LIST;
    }
    
    lapack_complex_double *result = (lapack_complex_double*)(stackpos+alen);
    LAPACK_zgeev("N", "N", &n, (lapack_complex_double*)stackpos, &lda, result, NULL, &ldv, NULL, &ldv, result+2*n, &lwork, (double*)(result+n), &info);
    lwork = creal(result[2*n]);
    LAPACK_zgeev("N", "N", &n, (lapack_complex_double*)stackpos, &lda, result, NULL, &ldv, NULL, &ldv, result+2*n, &lwork, (double*)(result+n), &info);
    for (int i=0; i < 2*n; i++) stackpos[i] = ((double*)result)[i]*SIGN_BIT(fs);
    return (n<<9) | TYPE_POINT;
}

uint32_t func_det(void *f, double *stackpos) {
    function *fs = (function*)f;
    function *arg = fs->first_arg;

    uint32_t argtype, alen, blen;
    argtype = arg->oper(arg, stackpos);
    alen = argtype>>8;
    int info;
    if (IS_TYPE(argtype, TYPE_POINT)) {
        // Complex matrix
        int32_t n = isqrt(alen/2);
        int *ipiv = (int*)(stackpos+alen);
        if (2*n*n != alen) FAIL("ERROR: can only compute determinant of square matrix\n");
        lapack_complex_double *a = (lapack_complex_double*)stackpos;
        LAPACK_zgetrf(&n, &n, a, &n, ipiv, &info);
        lapack_complex_double det = 1;
        for (int i=0; i < n; i++) {
            det *= a[(n+1)*i];
            if (i+1 != ipiv[i]) det = -det;
        }
        stackpos[0] = creal(det);
        stackpos[1] = cimag(det);
        return (1<<9) | TYPE_POINT;
    } else {
        int32_t n = isqrt(alen);
        int *ipiv = (int*)(stackpos+alen);
        if (n*n != alen) FAIL("ERROR: can only compute determinant of square matrix\n");
        LAPACK_dgetrf(&n, &n, stackpos, &n, ipiv, &info);
        double det = 1;
        for (int i=0; i < n; i++) {
            det *= stackpos[(n+1)*i];
            if (i+1 != ipiv[i]) det = -det;
        }
        stackpos[0] = det;
        return 1<<8;
    }
}

