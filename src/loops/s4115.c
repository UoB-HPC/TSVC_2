#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.11

//int s4115(int* __restrict__ ip)
real_t s4115(struct args_t * func_args)
{

//    indirect addressing
//    sparse dot product
//    gather is required

    int * __restrict__ ip = func_args->arg_info;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    real_t sum;
    for (int nl = 0; nl < iterations; nl++) {
        sum = 0.;
        for (int i = 0; i < LEN_1D; i++) {
            sum += a[i] * b[ip[i]];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return sum;
}

test_function_t p_s4115 = &s4115;
