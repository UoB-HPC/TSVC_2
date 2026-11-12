#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %5.1

real_t vif(struct args_t * func_args)
{

//    control loops
//    vector if

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            if (b[i] > (real_t)0.) {
                a[i] = b[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_vif = &vif;
