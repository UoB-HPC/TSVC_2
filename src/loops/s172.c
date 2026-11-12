#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %1.7

//int s172( int n1, int n3)
real_t s172(struct args_t * func_args)
{
//    symbolics
//    vectorizable if n3 .ne. 0

    struct{int a;int b;} * x = func_args->arg_info;
    int n1 = x->a;
    int n3 = x->b;

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int nl = 0; nl < iterations; nl++) {
        for (int i = n1-1; i < LEN_1D; i += n3) {
            a[i] += b[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
    return calc_checksum(__func__);
}

test_function_t p_s172 = &s172;
