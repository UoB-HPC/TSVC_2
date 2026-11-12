#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

real_t s2233(struct args_t * func_args)
{

//    loop interchange
//    interchanging with one of two inner loops

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 100*(iterations/LEN_2D); nl++) {
        for (int i = 1; i < LEN_2D; i++) {
            for (int j = 1; j < LEN_2D; j++) {
                aa[j][i] = aa[j-1][i] + cc[j][i];
            }
            for (int j = 1; j < LEN_2D; j++) {
                bb[i][j] = bb[i-1][j] + cc[i][j];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s2233 = &s2233;
