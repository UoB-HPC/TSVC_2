#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %4.12

real_t f(real_t a, real_t b){
    return a*b;
}

real_t s4121(struct args_t * func_args)
{

//    statement functions
//    elementwise multiplication

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            a[i] += f(b[i],c[i]);
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s4121 = &s4121;
