#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

real_t s1112(struct args_t * func_args)
{

//    linear dependence testing
//    loop reversal

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations*3; nl++) {
        for (int i = LEN_1D - 1; i >= 0; i--) {
            a[i] = b[i] + (real_t) 1.;
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s1112 = &s1112;

