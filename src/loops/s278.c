#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %2.7

real_t s278(struct args_t * func_args)
{

//    control flow
//    if/goto to block if-then-else

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            if (a[i] > (real_t)0.) {
                goto L20;
            }
            b[i] = -b[i] + d[i] * e[i];
            goto L30;
L20:
            c[i] = -c[i] + d[i] * e[i];
L30:
            a[i] = b[i] + c[i] * d[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s278 = &s278;
