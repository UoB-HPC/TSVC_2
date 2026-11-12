#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// **********************************************************
// *                                *
// *            VECTORIZATION                *
// *                                *
// **********************************************************

// %2.1

real_t s211(struct args_t * func_args)
{

//    statement reordering
//    statement reordering allows vectorization

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations; nl++) {
        for (int i = 1; i < LEN_1D-1; i++) {
            a[i] = b[i - 1] + c[i] * d[i];
            b[i] = b[i + 1] - e[i] * d[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s211 = &s211;
