#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %1.6

//int s162(int k)
real_t s162(struct args_t * func_args)
{
//    control flow
//    deriving assertions

    int k = *(int*)func_args->arg_info;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations; nl++) {
        if (k > 0) {
            for (int i = 0; i < LEN_1D-1; i++) {
                a[i] = a[i + k] + b[i] * c[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s162 = &s162;
