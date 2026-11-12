#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %1.3

real_t s131(struct args_t * func_args)
{
//    global data flow analysis
//    forward substitution

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    int m  = 1;
    for (int nl = 0; nl < 5*iterations; nl++) {
        for (int i = 0; i < LEN_1D - 1; i++) {
            a[i] = a[i + m] + b[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s131 = &s131;
