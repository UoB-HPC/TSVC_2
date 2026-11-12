#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

real_t s1421(struct args_t * func_args)
{

//    storage classes and equivalencing
//    equivalence- no overlap

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    xx = &b[LEN_1D/2];

    for (int nl = 0; nl < 8*iterations; nl++) {
        for (int i = 0; i < LEN_1D/2; i++) {
            b[i] = xx[i] + a[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 1.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s1421 = &s1421;
