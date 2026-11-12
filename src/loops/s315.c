#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// %3.1

real_t s315(struct args_t * func_args)
{

//    reductions
//    if to max with index reductio 1 dimension

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    for (int i = 0; i < LEN_1D; i++)
        a[i] = (i * 7) % LEN_1D;

    real_t x, chksum;
    int index;
    for (int nl = 0; nl < iterations; nl++) {
        x = a[0];
        index = 0;
        for (int i = 0; i < LEN_1D; ++i) {
            if (a[i] > x) {
                x = a[i];
                index = i;
            }
        }
        chksum = x + (real_t) index;
        dummy(a, b, c, d, e, aa, bb, cc, chksum);
    }

    gettimeofday(&func_args->t2, NULL);
    return index + x + 1;
}

test_function_t p_s315 = &s315;
