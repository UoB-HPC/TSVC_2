#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %2.5

real_t s254(struct args_t * func_args)
{

//    scalar and array expansion
//    carry around variable

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    real_t x;
    for (int nl = 0; nl < 4*iterations; nl++) {
        x = b[LEN_1D-1];
        for (int i = 0; i < LEN_1D; i++) {
            a[i] = (b[i] + x) * (real_t).5;
            x = b[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s254 = &s254;
