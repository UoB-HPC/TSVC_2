#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %1.1

real_t s113(struct args_t * func_args)
{

//    linear dependence testing
//    a(i)=a(1) but no actual dependence cycle

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 4*iterations; nl++) {
        for (int i = 1; i < LEN_1D; i++) {
            a[i] = a[0] + b[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s113 = &s113;

