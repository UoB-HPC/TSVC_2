#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %2.7

real_t s2712(struct args_t * func_args)
{

//    control flow
//    if to elemental min

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 4*iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            if (a[i] > b[i]) {
                a[i] += b[i] * c[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s2712 = &s2712;
