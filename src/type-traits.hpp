#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <concepts>
#include <random>
#include <array>

template <std::size_t size>
struct int_least {
    static_assert(size <= 128, "size must be less than or equal to 128");
    using type = typename std::conditional<
        size <= 8, std::int_least8_t,
        typename std::conditional<
            size <= 16, std::int_least16_t,
            typename std::conditional<
                size <= 32, std::int_least32_t,
                typename std::conditional<size <= 64, std::int_least64_t, __int128_t>::type>::type>::type>::type;
};

template <std::size_t size>
using int_least_t = typename int_least<size>::type;

template <std::size_t size>
struct uint_least {
    static_assert(size <= 128, "size must be less than or equal to 128");
    using type = typename std::conditional<
        size <= 8, std::uint_least8_t,
        typename std::conditional<
            size <= 16, std::uint_least16_t,
            typename std::conditional<
                size <= 32, std::uint_least32_t,
                typename std::conditional<size <= 64, std::uint_least64_t, __uint128_t>::type>::type>::type>::type;
};

template <typename T>
using signed_type_t = int_least_t<std::numeric_limits<T>::digits>;

using Torus16 = uint16_t;
using Torus = uint32_t;
Torus16 torus_to_torus16(Torus a) { return static_cast<Torus16>(a >> 16); }

template <typename G>
concept RandGen = std::uniform_random_bit_generator<G>;
