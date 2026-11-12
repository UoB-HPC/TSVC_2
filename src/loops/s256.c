#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %2.5

real_t s256(struct args_t * func_args)
{

//    scalar and array expansion
//    array expansion

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 10*(iterations/LEN_2D); nl++) {
        for (int i = 0; i < LEN_2D; i++) {
            for (int j = 1; j < LEN_2D; j++) {
                a[j] = (real_t)1.0 - a[j - 1];
                aa[j][i] = a[j] + bb[j][i]*d[j];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s256 = &s256;
