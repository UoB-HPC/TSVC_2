#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %1.5

void s151s(real_t a[LEN_1D], real_t b[LEN_1D],  int m)
{
    for (int i = 0; i < LEN_1D-1; i++) {
        a[i] = a[i + m] + b[i];
    }
}

real_t s151(struct args_t * func_args)
{

//    interprocedural data flow analysis
//    passing parameter information into a subroutine

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < 5*iterations; nl++) {
        s151s(a, b,  1);
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s151 = &s151;
