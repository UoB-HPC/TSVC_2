#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %2.9

real_t s291(struct args_t * func_args)
{

//    loop peeling
//    wrap around variable, 1 level

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    int im1;
    for (int nl = 0; nl < 2*iterations; nl++) {
        im1 = LEN_1D-1;
        for (int i = 0; i < LEN_1D; i++) {
            a[i] = (b[i] + b[im1]) * (real_t).5;
            im1 = i;
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s291 = &s291;
