#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %3.1

real_t s314(struct args_t * func_args)
{

//    reductions
//    if to max reduction

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    real_t x;
    for (int nl = 0; nl < iterations*5; nl++) {
        x = a[0];
        for (int i = 0; i < LEN_1D; i++) {
            if (a[i] > x) {
                x = a[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, x);
    }

    gettimeofday(&func_args->t2, NULL);
    return x;
}

test_function_t p_s314 = &s314;
