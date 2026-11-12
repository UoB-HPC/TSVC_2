#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.7

int s471s(void)
{
// --  dummy subroutine call made in s471
    return 0;
}

real_t s471(struct args_t * func_args){

//    call statements

    int m = LEN_1D;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations/2; nl++) {
        for (int i = 0; i < m; i++) {
            x[i] = b[i] + d[i] * d[i];
            s471s();
            b[i] = c[i] + d[i] * e[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s471 = &s471;
