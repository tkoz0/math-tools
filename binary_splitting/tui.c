// tkoz unsigned int

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// digits for conversion to string
const char *DIGITS = "0123456789abcdefghijklmnopqrstuvwxyz";

// unsigned integer data type
typedef struct
{
    uint32_t *arr;
    size_t len, alloc;
}
tui_t;

// 64 bit multiply to 128 bit result
static inline uint64_t mul64hi(uint64_t a, uint64_t b)
{
    return (a * (__uint128_t) b) >> 64;
}
static inline uint64_t mul64lo(uint64_t a, uint64_t b)
{
    return a * b;
}

// initialize integer
void tui_init(tui_t *a)
{
    a->arr = malloc(sizeof(*(a->arr)));
    assert(a->arr);
    a->len = 0;
    a->alloc = 1;
}

// set from 32 bit integer
void tui_set_ui(tui_t *a, uint32_t n)
{
    if (n == 0)
        a->len = 0;
    else
    {
        a->len = 1;
        a->arr[0] = n;
    }
}

// deallocate
void tui_free(tui_t *a)
{
    free(a->arr);
}

// change allocated space
void tui_resize(tui_t *a, size_t len)
{
    uint32_t *newptr = realloc(a->arr,len*sizeof(*(a->arr)));
    assert(newptr);
    for (size_t i = a->len; i < len; ++i)
        newptr[i] = 0;
    a->arr = newptr;
    a->alloc = len;
    if (len < a->len) // reduced size
    {
        a->len = len;
        while (a->len && a->arr[a->len-1] == 0)
            --(a->len);
    }
}

// set value from base 2^32 array
void tui_set_arr(tui_t *a, uint32_t *arr, size_t len)
{
    size_t nlen = 0;
    for (size_t i = 0; i < len; ++i)
        if (arr[i] != 0)
            nlen = i+1;
    if (a->alloc < nlen)
        tui_resize(a,nlen);
    memcpy(a->arr,arr,nlen*sizeof(*arr));
    a->len = nlen;
}

// TODO optimize for power of 2 base
char *tui_print(const tui_t *a, uint32_t base)
{
    assert(base >= 2 && base <= 36);
    size_t digits = (size_t)floor((a->len * 32) * log(2) / log(base)) + 1;
    char *ret = malloc(digits+1);
    if (a->len == 0) // special case for zero
    {
        memcpy(ret,"0",2);
        return ret;
    }
    uint32_t *buf = malloc(a->len*sizeof(*buf));
    size_t buflen = a->len;
    memcpy(buf,a->arr,a->len*sizeof(*buf));
    char *ptr = ret;
    while (buflen) // repeatedly divide by base
    {
        uint64_t value = 0;
        for (size_t i = buflen; i--;)
        {
            value <<= 32;
            value += buf[i];
            buf[i] = value / base;
            value %= base;
        }
        while (buflen && buf[buflen-1] == 0)
            --buflen;
        *(ptr++) = DIGITS[value]; // use remainder
    }
    free(buf);
    *ptr = '\0'; // terminate
    // reverse string order
    --ptr;
    for (char *ptr2 = ret; ptr2 < ptr; ++ptr2, --ptr)
    {
        char tmp = *ptr;
        *ptr = *ptr2;
        *ptr2 = tmp;
    }
    return ret;
}

// add 2 numbers
tui_t tui_add(const tui_t *restrict a, const tui_t *restrict b)
{
    size_t req_len = (a->len < b->len ? b->len : a->len) + 1; // max
    tui_t ret;
    tui_init(&ret);
    tui_resize(&ret,req_len);
    bool carry = 0;
    size_t i = 0, minlen = a->len < b->len ? a->len : b->len;
    while (i < minlen)
    {
        uint32_t sum = a->arr[i] + b->arr[i] + carry;
        carry = (sum < a->arr[i]) || (carry && sum == a->arr[i]);
        ret.arr[i++] = sum;
    }
    while (i < a->len) // use 0 in b
    {
        uint32_t sum = a->arr[i] + carry;
        carry = sum == 0;
        ret.arr[i++] = sum;
    }
    while (i < b->len) // use 0 in a
    {
        uint32_t sum = b->arr[i] + carry;
        carry = sum == 0;
        ret.arr[i++] = sum;
    }
    if (carry)
        ret.arr[i++] = 1;
    ret.len = i;
    return ret;
}

// multiply 2 numbers with the grid method and diagonal sums
// theoretical limit to max(a->len,b->len) < 2^32-1
// due to the use of 32 bit integer for higher bits of diagonal sums
tui_t tui_mul_basic_grid(const tui_t *restrict a, const tui_t *restrict b)
{
    size_t req_len = a->len + b->len;
    tui_t ret;
    tui_init(&ret);
    tui_resize(&ret,req_len);
    // compute the diagonal sums, using extra 32 bits for overflow
    uint64_t *sum0 = calloc(req_len-1,sizeof(*sum0));
    uint32_t *sum1 = calloc(req_len-1,sizeof(*sum1));
    // iterate over the full multiplication grid without storing it
    for (size_t i = 0; i < a->len; ++i)
        for (size_t j = 0; j < b->len; ++j)
        {
            uint64_t grid = (uint64_t)(a->arr[i]) * (uint64_t)(b->arr[j]);
            uint64_t sum = sum0[i+j] + grid;
            if (sum < grid)
                ++sum1[i+j];
            sum0[i+j] = sum;
        }
    uint64_t s0 = 0;
    uint32_t s1 = 0;
    for (size_t i = 0; i < req_len-1; ++i)
    {
        // add next element of diagonal sums to current sum
        uint64_t s0tmp = s0 + sum0[i];
        if (s0tmp < s0)
            ++s1;
        s1 += sum1[i]; // assumes no overflow
        // store lowest 32 bits and shift for next iteration
        ret.arr[i] = s0tmp;
        s0 = ((uint64_t)(s1) << 32) + (s0tmp >> 32);
        s1 = 0;
    }
    ret.arr[req_len-1] = s0;
    assert((s0 >> 32) == 0);
    ret.len = req_len;
    while (ret.len && ret.arr[ret.len-1] == 0)
        --(ret.len);
    free(sum0);
    free(sum1);
    return ret;
}

// multiply 2 numbers with the "grade school" algorithm 1 row at a time
tui_t tui_mul_basic_1row(const tui_t *restrict a, const tui_t *restrict b)
{
    size_t req_len = a->len + b->len;
    if (a->len < b->len) // ensure a is the longer one
    {
        const tui_t *tmp = a;
        a = b;
        b = tmp;
    }
    tui_t ret;
    tui_init(&ret);
    tui_resize(&ret,req_len);
    uint32_t *buf = malloc((b->len+1)*sizeof(*buf));
    for (size_t i = 0; i < a->len; ++i)
    {
        // compute a row of the "grade school method"
        uint32_t mul = a->arr[i];
        uint32_t mcarry = 0;
        for (size_t j = 0; j < b->len; ++j)
        {
            uint64_t prod = (uint64_t)(mul) * (uint64_t)(b->arr[j]) + mcarry;
            buf[j] = prod;
            mcarry = prod >> 32;
        }
        buf[b->len] = mcarry;
        // add row to the total
        bool acarry = 0;
        size_t j = i;
        for (; j < i + b->len + 1; ++j)
        {
            uint32_t sum = ret.arr[j] + buf[j-i] + acarry;
            acarry = (sum < ret.arr[j]) || (acarry && sum == ret.arr[j]);
            ret.arr[j] = sum;
        }
        while (acarry) // complete addition
        {
            ++ret.arr[j];
            acarry = ret.arr[j++] == 0;
        }
    }
    ret.len = req_len;
    while (ret.len && ret.arr[ret.len-1] == 0)
        --(ret.len);
    free(buf);
    return ret;
}

int main(int argc, char **argv)
{
    uint32_t a[3] = {57185,671,300};
    uint32_t b[3] = {91857201,5862865,999999999};
    tui_t ia, ib;
    tui_init(&ia);
    tui_init(&ib);
    tui_set_arr(&ia,a,3);
    tui_set_arr(&ib,b,3);
    tui_t ic = tui_mul_basic_1row(&ia,&ib);
    tui_t id = tui_add(&ia,&ib);
    char *sa = tui_print(&ia,10);
    char *sb = tui_print(&ib,10);
    char *sc = tui_print(&ic,10);
    char *sd = tui_print(&id,10);
    printf("%s\n%s\n%s\n%s\n",sa,sb,sc,sd);
    free(sa);
    free(sb);
    free(sc);
    free(sd);
    tui_free(&ia);
    tui_free(&ib);
    tui_free(&ic);
    tui_free(&id);
    return 0;
}
