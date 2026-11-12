#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

real_t s1251(struct args_t * func_args)
{

//    scalar and array expansion
//    scalar expansion

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    real_t s;
    for (int nl = 0; nl < 4*iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            s = b[i]+c[i];
            b[i] = a[i]+d[i];
            a[i] = s*e[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s1251 = &s1251;
