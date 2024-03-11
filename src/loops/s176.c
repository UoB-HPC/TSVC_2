#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %1.7

real_t s176(struct args_t * func_args)
{

//    symbolics
//    convolution

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    int m = LEN_1D/2;
    for (int nl = 0; nl < 4*(iterations/LEN_1D); nl++) {
        for (int j = 0; j < (LEN_1D/2); j++) {
            for (int i = 0; i < m; i++) {
                a[i] += b[i+m-j-1] * c[j];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s176 = &s176;
