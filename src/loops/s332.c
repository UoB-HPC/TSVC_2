#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %3.3
//int s332( real_t t)
real_t s332(struct args_t * func_args)
{

//    search loops
//    first value greater than threshold

    int t = *(int*)func_args->arg_info;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    int index;
    real_t value;
    real_t chksum;
    for (int nl = 0; nl < iterations; nl++) {
        index = -2;
        value = -1.;
        for (int i = 0; i < LEN_1D; i++) {
            if (a[i] > t) {
                index = i;
                value = a[i];
                goto L20;
            }
        }
L20:
        chksum = value + (real_t) index;
        dummy(a, b, c, d, e, aa, bb, cc, chksum);
    }

    gettimeofday(&func_args->t2, NULL);
    return value;
}

test_function_t p_s332 = &s332;
