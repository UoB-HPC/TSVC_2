#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// **********************************************************
//                                *
//             LANGUAGE COMPLETENESS            *
//                                *
// **********************************************************

// %4.1
// %4.2

real_t s421(struct args_t * func_args)
{

//    storage classes and equivalencing
//    equivalence- no overlap

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    xx = flat_2d_array;

    for (int nl = 0; nl < 4*iterations; nl++) {
        yy = xx;
        for (int i = 0; i < LEN_1D - 1; i++) {
            xx[i] = yy[i+1] + a[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 1.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s421 = &s421;
