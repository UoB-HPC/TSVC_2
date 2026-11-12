#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %2.7

//int s272(real_t t)
real_t s272(struct args_t * func_args)
{

//    control flow
//    loop with independent conditional

    int t = *(int*)func_args->arg_info;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            if (e[i] >= t) {
                a[i] += c[i] * d[i];
                b[i] += c[i] * c[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s272 = &s272;
