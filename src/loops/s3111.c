#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %3.1

real_t s3111(struct args_t * func_args)
{

//    reductions
//    conditional sum reduction

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    real_t sum;
    for (int nl = 0; nl < iterations/2; nl++) {
        sum = 0.;
        for (int i = 0; i < LEN_1D; i++) {
            if (a[i] > (real_t)0.) {
                sum += a[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, sum);
    }

    gettimeofday(&func_args->t2, NULL);
    return sum;
}

test_function_t p_s3111 = &s3111;
