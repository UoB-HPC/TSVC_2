#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.11

//int s4114(int* ip, int n1)
real_t s4114(struct args_t * func_args)
{

//    indirect addressing
//    mix indirect addressing with variable lower and upper bounds
//    gather is required

    struct{int * __restrict__ a;int b;} * x = func_args->arg_info;
    int * __restrict__ ip = x->a;
    int n1 = x->b;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    int k;
    for (int nl = 0; nl < iterations; nl++) {
        for (int i = n1-1; i < LEN_1D; i++) {
            k = ip[i];
            a[i] = b[i] + c[LEN_1D-k+1-2] * d[i];
            k += 5;
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s4114 = &s4114;
