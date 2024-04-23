#pragma once

#include <cstdint>

// compile time integer exponentiation

template <uint64_t b, uint64_t p> struct template_pow
{
    static const uint64_t value = b * template_pow<b,p-1>::value;
};
template <uint64_t b> struct template_pow<b,0>
{
    static const uint64_t value = 1;
};

#define TEMPLATE_POW(b,p) (template_pow<(b),(p)>::value)

template <uint64_t a, uint64_t b> struct template_gcd
{
    static_assert(a > 0 and b > 0);
    static const uint64_t value = template_gcd<b,a%b>::value;
};
template <uint64_t a> struct template_gcd<a,0>
{
    static_assert(a > 0);
    static const uint64_t value = a;
};
template <uint64_t a> struct template_gcd<0,a>
{
    static_assert(a > 0);
    static const uint64_t value = a;
};

#define TEMPLATE_GCD(a,b) (template_gcd<(a),(b)>::value)
