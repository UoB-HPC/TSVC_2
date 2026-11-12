#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.9

//int s491(int* __restrict__ ip)
real_t s491(struct args_t * func_args)
{

//    vector semantics
//    indirect addressing on lhs, store in sequence
//    scatter is required

    int * __restrict__ ip = func_args->arg_info;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            a[ip[i]] = b[i] + c[i] * d[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s491 = &s491;
