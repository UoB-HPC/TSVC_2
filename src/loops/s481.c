#include <stddef.h>
#include <stdlib.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.8

real_t s481(struct args_t * func_args)
{

//    non-local goto's
//    stop statement

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            if (d[i] < (real_t)0.) {
                exit (0);
            }
            a[i] += b[i] * c[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s481 = &s481;
