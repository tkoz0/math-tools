#include "bu64.hpp"

// high bits of 64 bit multiplication to 128 bit result
static inline uint64_t _mul64hi(uint64_t a, uint64_t b)
{ return ((__uint128_t) a * (__uint128_t) b) >> 64; }

// low bits of 64 bit multiplication to 128 bit result
static inline uint64_t _mul64lo(uint64_t a, uint64_t b)
{ return a * b; }

// full 128 bit result of 64 bit multiplication
typedef struct { uint64_t lo, hi; } _mul64_t;

// should compile to 3 instructions with -O3
static inline _mul64_t _mul64(uint64_t a, uint64_t b)
{
    _mul64_t ret;
    __uint128_t c = (__uint128_t) a * (__uint128_t) b;
    ret.lo = c;
    ret.hi = c >> 64;
    return ret;
}

// multiply big uint by uint64_t
uint64_t bu64_mul64(uint64_t *input, size_t ilen, uint64_t m)
{
    _mul64_t mv;
    uint64_t c = 0; // carry
    size_t i;
    for (i = 0; i < ilen; ++i)
    {
        mv = _mul64(m,input[i]);
        uint64_t tmp = mv.lo + c;
        input[i] = tmp;
        c = mv.hi + (tmp < mv.lo); // handle overflow
    }
    return c;
}

// add uint64_t to big uint
bool bu64_add64(uint64_t *input, size_t ilen, uint64_t a)
{
    uint64_t c = a;
    size_t i = 0;
    while (c and i < ilen)
    {
        uint64_t tmp = input[i] + c;
        input[i++] = tmp;
        c = (tmp < c);
    }
    return (i == ilen) and c;
}

// divide big uint by uint32_t (uint64_t does not k as easily)
uint32_t bu64_div32(uint64_t *input, size_t ilen, uint32_t d)
{
    uint32_t *input32 = (uint32_t*)input;
    uint64_t v = 0;
    for (size_t i = 2*ilen; i--;)
    {
        v = (v << 32) | input32[i]; // < d * 2^32
        input32[i] = v / d;
        v %= d;
    }
    return v;
}

// subtract uint64_t from big uint
bool bu64_sub64(uint64_t *input, size_t ilen, uint64_t a)
{
    uint64_t s = a;
    size_t i = 0;
    while (s and i < ilen)
    {
        uint64_t tmp = input[i];
        input[i] -= s;
        s = (input[i++] > tmp);
    }
    return not s;
}

const char *_digits = "0123456789abcdefghijklmnopqrstuvwxyz";

// character to numeric value for base conversion
static inline uint8_t _conv_digit(char c)
{
    return c <= '9' ? c-'0' : c-'a'+10;
}

// convert binary big uint to string
size_t bu64_write_str(uint8_t base, uint64_t *input, size_t ilen, char *output)
{
    char *optr = output;
    while (ilen and input[ilen-1] == 0)
        --ilen;
    while (ilen)
    {
        *(optr++) = _digits[bu64_div32(input,ilen,base)];
        if (input[ilen-1] == 0)
            --ilen;
    }
    *optr = '\0';
    size_t ret = optr - output;
    --optr;
    while (output < optr)
    {
        char tmp = *output;
        *(output++) = *optr;
        *(optr--) = tmp;
    }
    return ret;
}

// convert string to binary big uint
size_t bu64_read_str(uint8_t base, const char *input, uint64_t *output)
{
    size_t olen = 0;
    while (*input)
    {
        uint8_t d = _conv_digit(*(input++));
        uint64_t cm = bu64_mul64(output,olen,base);
        if (cm)
            output[olen++] = cm;
        bool ca = bu64_add64(output,olen,d);
        if (ca)
            output[olen++] = 1;
    }
    return olen;
}

#if 1 // === TESTING CODE BEGIN ===

#include <cstdio>
#include <vector>
#define base 10
int main(int argc, char **argv)
{
    std::vector<uint64_t> arr;
    arr.push_back(1);
    for (uint64_t i = 1; i <= (uint64_t)atoll(argv[1]); ++i)
    {
        uint64_t tmp = bu64_mul64(arr.data(),arr.size(),i);
        if (tmp)
            arr.push_back(tmp);
    }
    size_t slen = 1+(size_t)((arr.size()<<6)/log2(base));
    printf("est str length = %lu\n",slen);
    char *b = (char*)malloc(slen);
    bu64_write_str(base,arr.data(),arr.size(),b);
    printf("out str length = %lu\n",strlen(b));
    //printf("%s\n",b);
    free(b);
    return 0;
};

#endif // === TESTING CODE END ===

/* possibly useful intrinsics (all values in "little endian" order)

[32 bit mult -> 64 bit result]
(AVX2) _mm256_mul_epu32([a0,0,a1,0,a2,0,a3,0],[b0,0,b1,0,b2,0,b3,0])
-> [a0*b0,a1*b1,a2*b2,a3*b3]
(SSE2) _mm_mul_epu32([a0,0,a1,0],[b0,0,b1,0]) -> [a0*b0,a1*b1]

[parallel addition/subtraction]
(AVX2) _mm256_sub_epi32([a0..a7],[b0..b7]) -> [ai-bi, ...]
(AVX2) _mm256_add_epi32([a0..a7],[b0..b7]) -> [ai+bi, ...]
(AVX2) _mm256_sub_epi64([a0..a3],[b0..b3]) -> [ai-bi, ...]
(AVX2) _mm256_add_epi64([a0..a3],[b0..b3]) -> [ai+bi, ...]
(SSE2) _mm_add_epi32([a0..a3],[b0..b3]) -> [ai+bi, ...]
(SSE2) _mm_add_epi64([a0,a1],[b0,b1]) -> [a0+b0,a1+b1]
(SSE2) _mm_sub_epi32([a0..a3],[b0..b3]) -> [ai-bi, ...]
(SSE2) _mm_sub_epi64([a0,a1],[b0,b1]) -> [a0-b0,a1-b1]

[set values]
(AVX) _mm256_set_epi32(a0,a1,a2,a3,a4,a5,a6,a7)
(AVX) _mm256_set_epi64x(a0,a1,a2,a3)
(SSE2) _mm_set_epi32(a0,a1,a2,a3)
(SSE2) _mm_set_epi64x(a0,a1)

*/
