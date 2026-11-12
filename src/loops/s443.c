#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.4

real_t s443(struct args_t * func_args)
{

//    non-logical if's
//    arithmetic if

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 2*iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            if (d[i] <= (real_t)0.) {
                goto L20;
            } else {
                goto L30;
            }
L20:
            a[i] += b[i] * c[i];
            goto L50;
L30:
            a[i] += b[i] * b[i];
L50:
            ;
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s443 = &s443;
