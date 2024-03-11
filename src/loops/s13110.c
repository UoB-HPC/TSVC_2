#include <stddef.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

real_t s13110(struct args_t * func_args)
{

//    reductions
//    if to max with index reductio 2 dimensions

    initialise_arrays(__func__);
    gettimeofday(&func_args->t1, NULL);

    int xindex, yindex;
    real_t max, chksum;
    for (int nl = 0; nl < 100*(iterations/(LEN_2D)); nl++) {
        max = aa[(0)][0];
        xindex = 0;
        yindex = 0;
        for (int i = 0; i < LEN_2D; i++) {
            for (int j = 0; j < LEN_2D; j++) {
                if (aa[i][j] > max) {
                    max = aa[i][j];
                    xindex = i;
                    yindex = j;
                }
            }
        }
        chksum = max + (real_t) xindex + (real_t) yindex;
        dummy(a, b, c, d, e, aa, bb, cc, chksum);
    }

    gettimeofday(&func_args->t2, NULL);
    return max + xindex+1 + yindex+1;
}

test_function_t p_s13110 = &s13110;
