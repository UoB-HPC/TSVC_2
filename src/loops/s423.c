#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.2

real_t s423(struct args_t * func_args)
{

//    storage classes and equivalencing
//    common and equivalenced variables - with anti-dependence

    // do this again here
    int vl = 64;
    xx = flat_2d_array + vl;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 4*iterations; nl++) {
        for (int i = 0; i < LEN_1D - 1; i++) {
            flat_2d_array[i+1] = xx[i] + a[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 1.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s423 = &s423;
