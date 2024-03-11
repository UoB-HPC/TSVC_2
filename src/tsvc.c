
/*
 * This is an executable test containing a number of loops to measure
 * the performance of a compiler. Arrays' length is LEN_1D by default
 * and if you want a different array length, you should replace every
 * LEN_1D by your desired number which must be a multiple of 40. If you
 * want to increase the number of loop calls to have a longer run time
 * you have to manipulate the constant value iterations. There is a dummy
 * function called in each loop to make all computations appear required.
 * The time to execute this function is included in the time measurement
 * for the output but it is neglectable.
 *
 *  The output includes three columns:
 *    Loop:        The name of the loop
 *    Time(Sec):     The time in seconds to run the loop
 *    Checksum:    The checksum calculated when the test has run
 *
 * In this version of the codelets arrays are static type.
 *
 * All functions/loops are taken from "TEST SUITE FOR VECTORIZING COMPILERS"
 * by David Callahan, Jack Dongarra and David Levine except those whose
 * functions' name have 4 digits.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"

// array definitions
__attribute__((aligned(ARRAY_ALIGNMENT))) real_t flat_2d_array[LEN_2D*LEN_2D];

__attribute__((aligned(ARRAY_ALIGNMENT))) real_t x[LEN_1D];

__attribute__((aligned(ARRAY_ALIGNMENT))) real_t a[LEN_1D],b[LEN_1D],c[LEN_1D],d[LEN_1D],e[LEN_1D],
                                   aa[LEN_2D][LEN_2D],bb[LEN_2D][LEN_2D],cc[LEN_2D][LEN_2D],tt[LEN_2D][LEN_2D];

__attribute__((aligned(ARRAY_ALIGNMENT))) int indx[LEN_1D];

real_t* __restrict__ xx;
real_t* yy;

// FIXME better spot for this (maybe macro-generate)
test_function_t p_s000 __attribute__((weak)) = NULL;
test_function_t p_s111 __attribute__((weak)) = NULL;
test_function_t p_s1111 __attribute__((weak)) = NULL;
test_function_t p_s112 __attribute__((weak)) = NULL;
test_function_t p_s1112 __attribute__((weak)) = NULL;
test_function_t p_s113 __attribute__((weak)) = NULL;
test_function_t p_s1113 __attribute__((weak)) = NULL;
test_function_t p_s114 __attribute__((weak)) = NULL;
test_function_t p_s115 __attribute__((weak)) = NULL;
test_function_t p_s1115 __attribute__((weak)) = NULL;
test_function_t p_s116 __attribute__((weak)) = NULL;
test_function_t p_s118 __attribute__((weak)) = NULL;
test_function_t p_s119 __attribute__((weak)) = NULL;
test_function_t p_s1119 __attribute__((weak)) = NULL;
test_function_t p_s121 __attribute__((weak)) = NULL;
test_function_t p_s122 __attribute__((weak)) = NULL;
test_function_t p_s123 __attribute__((weak)) = NULL;
test_function_t p_s124 __attribute__((weak)) = NULL;
test_function_t p_s125 __attribute__((weak)) = NULL;
test_function_t p_s126 __attribute__((weak)) = NULL;
test_function_t p_s127 __attribute__((weak)) = NULL;
test_function_t p_s128 __attribute__((weak)) = NULL;
test_function_t p_s131 __attribute__((weak)) = NULL;
test_function_t p_s132 __attribute__((weak)) = NULL;
test_function_t p_s141 __attribute__((weak)) = NULL;
test_function_t p_s151 __attribute__((weak)) = NULL;
test_function_t p_s152 __attribute__((weak)) = NULL;
test_function_t p_s161 __attribute__((weak)) = NULL;
test_function_t p_s1161 __attribute__((weak)) = NULL;
test_function_t p_s162 __attribute__((weak)) = NULL;
test_function_t p_s171 __attribute__((weak)) = NULL;
test_function_t p_s172 __attribute__((weak)) = NULL;
test_function_t p_s173 __attribute__((weak)) = NULL;
test_function_t p_s174 __attribute__((weak)) = NULL;
test_function_t p_s175 __attribute__((weak)) = NULL;
test_function_t p_s176 __attribute__((weak)) = NULL;
test_function_t p_s211 __attribute__((weak)) = NULL;
test_function_t p_s212 __attribute__((weak)) = NULL;
test_function_t p_s1213 __attribute__((weak)) = NULL;
test_function_t p_s221 __attribute__((weak)) = NULL;
test_function_t p_s1221 __attribute__((weak)) = NULL;
test_function_t p_s222 __attribute__((weak)) = NULL;
test_function_t p_s231 __attribute__((weak)) = NULL;
test_function_t p_s232 __attribute__((weak)) = NULL;
test_function_t p_s1232 __attribute__((weak)) = NULL;
test_function_t p_s233 __attribute__((weak)) = NULL;
test_function_t p_s2233 __attribute__((weak)) = NULL;
test_function_t p_s235 __attribute__((weak)) = NULL;
test_function_t p_s241 __attribute__((weak)) = NULL;
test_function_t p_s242 __attribute__((weak)) = NULL;
test_function_t p_s243 __attribute__((weak)) = NULL;
test_function_t p_s244 __attribute__((weak)) = NULL;
test_function_t p_s1244 __attribute__((weak)) = NULL;
test_function_t p_s2244 __attribute__((weak)) = NULL;
test_function_t p_s251 __attribute__((weak)) = NULL;
test_function_t p_s1251 __attribute__((weak)) = NULL;
test_function_t p_s2251 __attribute__((weak)) = NULL;
test_function_t p_s3251 __attribute__((weak)) = NULL;
test_function_t p_s252 __attribute__((weak)) = NULL;
test_function_t p_s253 __attribute__((weak)) = NULL;
test_function_t p_s254 __attribute__((weak)) = NULL;
test_function_t p_s255 __attribute__((weak)) = NULL;
test_function_t p_s256 __attribute__((weak)) = NULL;
test_function_t p_s257 __attribute__((weak)) = NULL;
test_function_t p_s258 __attribute__((weak)) = NULL;
test_function_t p_s261 __attribute__((weak)) = NULL;
test_function_t p_s271 __attribute__((weak)) = NULL;
test_function_t p_s272 __attribute__((weak)) = NULL;
test_function_t p_s273 __attribute__((weak)) = NULL;
test_function_t p_s274 __attribute__((weak)) = NULL;
test_function_t p_s275 __attribute__((weak)) = NULL;
test_function_t p_s2275 __attribute__((weak)) = NULL;
test_function_t p_s276 __attribute__((weak)) = NULL;
test_function_t p_s277 __attribute__((weak)) = NULL;
test_function_t p_s278 __attribute__((weak)) = NULL;
test_function_t p_s279 __attribute__((weak)) = NULL;
test_function_t p_s1279 __attribute__((weak)) = NULL;
test_function_t p_s2710 __attribute__((weak)) = NULL;
test_function_t p_s2711 __attribute__((weak)) = NULL;
test_function_t p_s2712 __attribute__((weak)) = NULL;
test_function_t p_s281 __attribute__((weak)) = NULL;
test_function_t p_s1281 __attribute__((weak)) = NULL;
test_function_t p_s291 __attribute__((weak)) = NULL;
test_function_t p_s292 __attribute__((weak)) = NULL;
test_function_t p_s293 __attribute__((weak)) = NULL;
test_function_t p_s2101 __attribute__((weak)) = NULL;
test_function_t p_s2102 __attribute__((weak)) = NULL;
test_function_t p_s2111 __attribute__((weak)) = NULL;
test_function_t p_s311 __attribute__((weak)) = NULL;
test_function_t p_s31111 __attribute__((weak)) = NULL;
test_function_t p_s312 __attribute__((weak)) = NULL;
test_function_t p_s313 __attribute__((weak)) = NULL;
test_function_t p_s314 __attribute__((weak)) = NULL;
test_function_t p_s315 __attribute__((weak)) = NULL;
test_function_t p_s316 __attribute__((weak)) = NULL;
test_function_t p_s317 __attribute__((weak)) = NULL;
test_function_t p_s318 __attribute__((weak)) = NULL;
test_function_t p_s319 __attribute__((weak)) = NULL;
test_function_t p_s3110 __attribute__((weak)) = NULL;
test_function_t p_s13110 __attribute__((weak)) = NULL;
test_function_t p_s3111 __attribute__((weak)) = NULL;
test_function_t p_s3112 __attribute__((weak)) = NULL;
test_function_t p_s3113 __attribute__((weak)) = NULL;
test_function_t p_s321 __attribute__((weak)) = NULL;
test_function_t p_s322 __attribute__((weak)) = NULL;
test_function_t p_s323 __attribute__((weak)) = NULL;
test_function_t p_s331 __attribute__((weak)) = NULL;
test_function_t p_s332 __attribute__((weak)) = NULL;
test_function_t p_s341 __attribute__((weak)) = NULL;
test_function_t p_s342 __attribute__((weak)) = NULL;
test_function_t p_s343 __attribute__((weak)) = NULL;
test_function_t p_s351 __attribute__((weak)) = NULL;
test_function_t p_s1351 __attribute__((weak)) = NULL;
test_function_t p_s352 __attribute__((weak)) = NULL;
test_function_t p_s353 __attribute__((weak)) = NULL;
test_function_t p_s421 __attribute__((weak)) = NULL;
test_function_t p_s1421 __attribute__((weak)) = NULL;
test_function_t p_s422 __attribute__((weak)) = NULL;
test_function_t p_s423 __attribute__((weak)) = NULL;
test_function_t p_s424 __attribute__((weak)) = NULL;
test_function_t p_s431 __attribute__((weak)) = NULL;
test_function_t p_s441 __attribute__((weak)) = NULL;
test_function_t p_s442 __attribute__((weak)) = NULL;
test_function_t p_s443 __attribute__((weak)) = NULL;
test_function_t p_s451 __attribute__((weak)) = NULL;
test_function_t p_s452 __attribute__((weak)) = NULL;
test_function_t p_s453 __attribute__((weak)) = NULL;
test_function_t p_s471 __attribute__((weak)) = NULL;
test_function_t p_s481 __attribute__((weak)) = NULL;
test_function_t p_s482 __attribute__((weak)) = NULL;
test_function_t p_s491 __attribute__((weak)) = NULL;
test_function_t p_s4112 __attribute__((weak)) = NULL;
test_function_t p_s4113 __attribute__((weak)) = NULL;
test_function_t p_s4114 __attribute__((weak)) = NULL;
test_function_t p_s4115 __attribute__((weak)) = NULL;
test_function_t p_s4116 __attribute__((weak)) = NULL;
test_function_t p_s4117 __attribute__((weak)) = NULL;
test_function_t p_s4121 __attribute__((weak)) = NULL;
test_function_t p_va __attribute__((weak)) = NULL;
test_function_t p_vag __attribute__((weak)) = NULL;
test_function_t p_vas __attribute__((weak)) = NULL;
test_function_t p_vif __attribute__((weak)) = NULL;
test_function_t p_vpv __attribute__((weak)) = NULL;
test_function_t p_vtv __attribute__((weak)) = NULL;
test_function_t p_vpvtv __attribute__((weak)) = NULL;
test_function_t p_vpvts __attribute__((weak)) = NULL;
test_function_t p_vpvpv __attribute__((weak)) = NULL;
test_function_t p_vtvtv __attribute__((weak)) = NULL;
test_function_t p_vsumr __attribute__((weak)) = NULL;
test_function_t p_vdotr __attribute__((weak)) = NULL;
test_function_t p_vbor __attribute__((weak)) = NULL;

void time_function(test_function_t vector_func, void * arg_info)
{
    if (!vector_func) return;

    struct args_t func_args = {.arg_info=arg_info};

    double result = vector_func(&func_args);

    double tic=func_args.t1.tv_sec+(func_args.t1.tv_usec/1000000.0);
    double toc=func_args.t2.tv_sec+(func_args.t2.tv_usec/1000000.0);

    double taken = toc-tic;

    printf("%10.3f\t%f\n", taken, result);
}

int main(void){
    int n1 = 1;
    int n3 = 1;
    int* ip;
    real_t s1,s2;
    init(&ip, &s1, &s2);
    printf("Loop \tTime(sec) \tChecksum\n");

    time_function(p_s000, NULL);
    time_function(p_s111, NULL);
    time_function(p_s1111, NULL);
    time_function(p_s112, NULL);
    time_function(p_s1112, NULL);
    time_function(p_s113, NULL);
    time_function(p_s1113, NULL);
    time_function(p_s114, NULL);
    time_function(p_s115, NULL);
    time_function(p_s1115, NULL);
    time_function(p_s116, NULL);
    time_function(p_s118, NULL);
    time_function(p_s119, NULL);
    time_function(p_s1119, NULL);
    time_function(p_s121, NULL);
    time_function(p_s122, &(struct{int a;int b;}){n1, n3});
    time_function(p_s123, NULL);
    time_function(p_s124, NULL);
    time_function(p_s125, NULL);
    time_function(p_s126, NULL);
    time_function(p_s127, NULL);
    time_function(p_s128, NULL);
    time_function(p_s131, NULL);
    time_function(p_s132, NULL);
    time_function(p_s141, NULL);
    time_function(p_s151, NULL);
    time_function(p_s152, NULL);
    time_function(p_s161, NULL);
    time_function(p_s1161, NULL);
    time_function(p_s162, &n1);
    time_function(p_s171, &n1);
    time_function(p_s172, &(struct{int a;int b;}){n1, n3});
    time_function(p_s173, NULL);
    time_function(p_s174, &(struct{int a;}){LEN_1D/2});
    time_function(p_s175, &n1);
    time_function(p_s176, NULL);
    time_function(p_s211, NULL);
    time_function(p_s212, NULL);
    time_function(p_s1213, NULL);
    time_function(p_s221, NULL);
    time_function(p_s1221, NULL);
    time_function(p_s222, NULL);
    time_function(p_s231, NULL);
    time_function(p_s232, NULL);
    time_function(p_s1232, NULL);
    time_function(p_s233, NULL);
    time_function(p_s2233, NULL);
    time_function(p_s235, NULL);
    time_function(p_s241, NULL);
    time_function(p_s242, &(struct{real_t a;real_t b;}){s1, s2});
    time_function(p_s243, NULL);
    time_function(p_s244, NULL);
    time_function(p_s1244, NULL);
    time_function(p_s2244, NULL);
    time_function(p_s251, NULL);
    time_function(p_s1251, NULL);
    time_function(p_s2251, NULL);
    time_function(p_s3251, NULL);
    time_function(p_s252, NULL);
    time_function(p_s253, NULL);
    time_function(p_s254, NULL);
    time_function(p_s255, NULL);
    time_function(p_s256, NULL);
    time_function(p_s257, NULL);
    time_function(p_s258, NULL);
    time_function(p_s261, NULL);
    time_function(p_s271, NULL);
    time_function(p_s272, &s1);
    time_function(p_s273, NULL);
    time_function(p_s274, NULL);
    time_function(p_s275, NULL);
    time_function(p_s2275, NULL);
    time_function(p_s276, NULL);
    time_function(p_s277, NULL);
    time_function(p_s278, NULL);
    time_function(p_s279, NULL);
    time_function(p_s1279, NULL);
    time_function(p_s2710, &s1);
    time_function(p_s2711, NULL);
    time_function(p_s2712, NULL);
    time_function(p_s281, NULL);
    time_function(p_s1281, NULL);
    time_function(p_s291, NULL);
    time_function(p_s292, NULL);
    time_function(p_s293, NULL);
    time_function(p_s2101, NULL);
    time_function(p_s2102, NULL);
    time_function(p_s2111, NULL);
    time_function(p_s311, NULL);
    time_function(p_s31111, NULL);
    time_function(p_s312, NULL);
    time_function(p_s313, NULL);
    time_function(p_s314, NULL);
    time_function(p_s315, NULL);
    time_function(p_s316, NULL);
    time_function(p_s317, NULL);
    time_function(p_s318, &n1);
    time_function(p_s319, NULL);
    time_function(p_s3110, NULL);
    time_function(p_s13110, NULL);
    time_function(p_s3111, NULL);
    time_function(p_s3112, NULL);
    time_function(p_s3113, NULL);
    time_function(p_s321, NULL);
    time_function(p_s322, NULL);
    time_function(p_s323, NULL);
    time_function(p_s331, NULL);
    time_function(p_s332, &s1);
    time_function(p_s341, NULL);
    time_function(p_s342, NULL);
    time_function(p_s343, NULL);
    time_function(p_s351, NULL);
    time_function(p_s1351, NULL);
    time_function(p_s352, NULL);
    time_function(p_s353, ip);
    time_function(p_s421, NULL);
    time_function(p_s1421, NULL);
    time_function(p_s422, NULL);
    time_function(p_s423, NULL);
    time_function(p_s424, NULL);
    time_function(p_s431, NULL);
    time_function(p_s441, NULL);
    time_function(p_s442, NULL);
    time_function(p_s443, NULL);
    time_function(p_s451, NULL);
    time_function(p_s452, NULL);
    time_function(p_s453, NULL);
    time_function(p_s471, NULL);
    time_function(p_s481, NULL);
    time_function(p_s482, NULL);
    time_function(p_s491, ip);
    time_function(p_s4112, &(struct{int*a;real_t b;}){ip, s1});
    time_function(p_s4113, ip);
    time_function(p_s4114, &(struct{int*a;int b;}){ip, n1});
    time_function(p_s4115, ip);
    time_function(p_s4116, &(struct{int * a; int b; int c;}){ip, LEN_2D/2, n1});
    time_function(p_s4117, NULL);
    time_function(p_s4121, NULL);
    time_function(p_va, NULL);
    time_function(p_vag, ip);
    time_function(p_vas, ip);
    time_function(p_vif, NULL);
    time_function(p_vpv, NULL);
    time_function(p_vtv, NULL);
    time_function(p_vpvtv, NULL);
    time_function(p_vpvts, &s1);
    time_function(p_vpvpv, NULL);
    time_function(p_vtvtv, NULL);
    time_function(p_vsumr, NULL);
    time_function(p_vdotr, NULL);
    time_function(p_vbor, NULL);

    return EXIT_SUCCESS;
}
