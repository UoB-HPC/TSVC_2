#include <stddef.h>
#include <math.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %3.1

real_t s3113(struct args_t * func_args)
{

//    reductions
//    maximum of absolute value

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    real_t max;
    for (int nl = 0; nl < iterations*4; nl++) {
        max = ABS(a[0]);
        for (int i = 0; i < LEN_1D; i++) {
            if ((ABS(a[i])) > max) {
                max = ABS(a[i]);
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, max);
    }

    gettimeofday(&func_args->t2, NULL);
    return max;
}

test_function_t p_s3113 = &s3113;
