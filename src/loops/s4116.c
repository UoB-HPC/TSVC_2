#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.11

//int s4116(int* __restrict__ ip, int j, int inc)
real_t s4116(struct args_t * func_args)
{

//    indirect addressing
//    more complicated sparse sdot
//    gather is required

    struct{int * __restrict__ a;int b;int c;} * x = func_args->arg_info;
    int * __restrict__ ip = x->a;
    int j = x->b;
    int inc = x->c;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    real_t sum;
    int off;
    for (int nl = 0; nl < 100*iterations; nl++) {
        sum = 0.;
        for (int i = 0; i < LEN_2D-1; i++) {
            off = inc + i;
            sum += a[off] * aa[j-1][ip[i]];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return sum;
}

test_function_t p_s4116 = &s4116;
