#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %3.1
real_t s313(struct args_t * func_args)
{

//    reductions
//    dot product

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    real_t dot;
    for (int nl = 0; nl < iterations*5; nl++) {
        dot = (real_t)0.;
        for (int i = 0; i < LEN_1D; i++) {
            dot += a[i] * b[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, dot);
    }

    gettimeofday(&func_args->t2, NULL);
    return dot;
}

test_function_t p_s313 = &s313;
