#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %2.4

//int s242(real_t s1, real_t s2)
real_t s242(struct args_t * func_args)
{

//    node splitting

    struct{real_t a;real_t b;} * x = func_args->arg_info;
    real_t s1 = x->a;
    real_t s2 = x->b;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations/5; nl++) {
        for (int i = 1; i < LEN_1D; ++i) {
            a[i] = a[i - 1] + s1 + s2 + b[i] + c[i] + d[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s242 = &s242;
