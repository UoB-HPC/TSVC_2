#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.5

real_t s453(struct args_t * func_args)
{

//    induction varibale recognition

    real_t s;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations*2; nl++) {
        s = 0.;
        for (int i = 0; i < LEN_1D; i++) {
            s += (real_t)2.;
            a[i] = s * b[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s453 = &s453;
