#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %1.1

real_t s112(struct args_t * func_args)
{

//    linear dependence testing
//    loop reversal

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 3*iterations; nl++) {
        for (int i = LEN_1D - 2; i >= 0; i--) {
            a[i+1] = a[i] + b[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s112 = &s112;

