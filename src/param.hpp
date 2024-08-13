#pragma once
#include <cstdint>

namespace param {
class Security128bit {
   public:
    static constexpr std::size_t n = 636;
    static constexpr double alpha = 0.0000925119974676756;
    static constexpr std::size_t Nbit = 9;
    static constexpr std::size_t N = 1 << Nbit;
    static constexpr std::size_t k = 2;
    static constexpr double alpha_bk = 0.0000000342348787018369;
    static constexpr std::size_t Bgbit = 8;
    static constexpr int Bg = 1 << Bgbit;
    static constexpr std::size_t l = 2;
};
}  // namespace param