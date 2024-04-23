// low level test implementation of large unsigned integer type
// representation using 64 bit limbs (uint64_t)
// function prefix bu64_ (big unsigned integer with 64 bit limbs)

#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>

// multiply by m, returns the last carry amount
// (nonzero last carry means the number gets longer)
uint64_t bu64_mul64(uint64_t *input, size_t ilen, uint64_t m);

// add a, return true if carry propagates past end
// (carry propagation past end means number gets longer)
bool bu64_add64(uint64_t *input, size_t ilen, uint64_t a);

// divide by d, return remainder
uint32_t bu64_div32(uint64_t *input, size_t ilen, uint32_t d);

// subtract a, return false if a is bigger that the input
bool bu64_sub64(uint64_t *input, size_t ilen, uint64_t a);

// convert number to string, output must be long enough (including null)
size_t bu64_write_str(uint8_t base, uint64_t *input, size_t ilen, char *output);

// convert string to number, output must be long enough
size_t bu64_read_str(uint8_t base, const char *input, uint64_t *output);

