#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.5

real_t s452(struct args_t * func_args)
{

//    intrinsic functions
//    seq function

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 4*iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            a[i] = b[i] + c[i] * (real_t) (i+1);
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s452 = &s452;
