#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %1.2

real_t s127(struct args_t * func_args)
{

//    induction variable recognition
//    induction variable with multiple increments

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    int j;
    for (int nl = 0; nl < 2*iterations; nl++) {
        j = -1;
        for (int i = 0; i < LEN_1D/2; i++) {
            j++;
            a[j] = b[i] + c[i] * d[i];
            j++;
            a[j] = b[i] + d[i] * e[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s127 = &s127;
