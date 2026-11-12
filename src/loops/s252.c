#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %2.5

real_t s252(struct args_t * func_args)
{

//    scalar and array expansion
//    loop with ambiguous scalar temporary

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    real_t t, s;
    for (int nl = 0; nl < iterations; nl++) {
        t = (real_t) 0.;
        for (int i = 0; i < LEN_1D; i++) {
            s = b[i] * c[i];
            a[i] = s + t;
            t = s;
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s252 = &s252;
