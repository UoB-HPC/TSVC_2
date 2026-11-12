#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

real_t s1111(struct args_t * func_args)
{

//    no dependence - vectorizable
//    jump in data access

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 2*iterations; nl++) {
        for (int i = 0; i < LEN_1D/2; i++) {
            a[2*i] = c[i] * b[i] + d[i] * b[i] + c[i] * c[i] + d[i] * b[i] + d[i] * c[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s1111 = &s1111;

