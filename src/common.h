
#define iterations 100000
#define LEN_1D 32000
#define LEN_2D 256

#define ARRAY_ALIGNMENT 64

#if 0
typedef double real_t;
#define ABS fabs
#else
typedef float real_t;
#define ABS fabsf
#endif

int dummy(real_t[LEN_1D], real_t[LEN_1D], real_t[LEN_1D], real_t[LEN_1D], real_t[LEN_1D], real_t[LEN_2D][LEN_2D], real_t[LEN_2D][LEN_2D], real_t[LEN_2D][LEN_2D], real_t);

void init(int** ip, real_t* s1, real_t* s2);

int initialise_arrays(const char* name);
real_t calc_checksum(const char * name);

