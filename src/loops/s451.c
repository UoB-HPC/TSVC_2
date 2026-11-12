#include <stddef.h>
#include <math.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.5

real_t s451(struct args_t * func_args)
{

//    intrinsic functions
//    intrinsics

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations/5; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            a[i] = sinf(b[i]) + cosf(c[i]);
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s451 = &s451;
