#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %5.1

//int vas( int* __restrict__ ip)
real_t vas(struct args_t * func_args)
{

//    control loops
//    vector assignment, scatter
//    scatter is required

    int * __restrict__ ip = func_args->arg_info;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 2*iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            a[ip[i]] = b[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_vas = &vas;
