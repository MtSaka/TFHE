#pragma once

#include "type-traits.hpp"
#include <random>

Torus16 double_to_torus16(double dbl) {
    return static_cast<Torus16>(std::fmod(dbl, 1.0) *
                                pow(2.0, std::numeric_limits<Torus16>::digits));
}

Torus double_to_torus(double dbl) {
    return static_cast<Torus>(std::fmod(dbl, 1.0) *
                              pow(2.0, std::numeric_limits<Torus>::digits));
}

template <RandGen Gen>
Torus16 uniform_torus16_gen(Gen& rng) {
    return std::uniform_int_distribution<Torus16>{
        0, std::numeric_limits<Torus16>::max()}(rng);
}

template <RandGen Gen>
Torus uniform_torus_gen(Gen& rng) {
    return std::uniform_int_distribution<Torus>{
        0, std::numeric_limits<Torus>::max()}(rng);
}

template <RandGen Gen>
Torus16 normal_torus16_gen(double alpha, Gen& rng) {
    return double_to_torus16(std::normal_distribution<>{0.0, alpha}(rng));
}
template <RandGen Gen>
Torus normal_torus_gen(double alpha, Gen& rng) {
    return double_to_torus(std::normal_distribution<>{0.0, alpha}(rng));
}

template <RandGen Gen>
uint32_t uniform_integer(uint32_t l, uint32_t r, Gen& rng) {
    return std::uniform_int_distribution<uint32_t>{l, r}(rng);
}