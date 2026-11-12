#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

real_t s1115(struct args_t * func_args)
{

//    linear dependence testing
//    triangular saxpy loop

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 100*(iterations/LEN_2D); nl++) {
        for (int i = 0; i < LEN_2D; i++) {
            for (int j = 0; j < LEN_2D; j++) {
                aa[i][j] = aa[i][j]*cc[j][i] + bb[i][j];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s1115 = &s1115;

