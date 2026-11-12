#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %1.1

real_t s119(struct args_t * func_args)
{

//    linear dependence testing
//    no dependence - vectorizable

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 200*(iterations/(LEN_2D)); nl++) {
        for (int i = 1; i < LEN_2D; i++) {
            for (int j = 1; j < LEN_2D; j++) {
                aa[i][j] = aa[i-1][j-1] + bb[i][j];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s119 = &s119;

