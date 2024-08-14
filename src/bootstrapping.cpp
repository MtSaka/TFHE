#pragma once

#include "trgsw.cpp"

#include <memory>

template <class Parameter>
struct BootstrappingKey {
   private:
    std::array<TRGSW<Parameter>, Parameter::n> data;

   public:
    BootstrappingKey() : data({}) {}
    BootstrappingKey(const BootstrappingKey& bk) : data(bk.data) {}
    template <RandGen Gen>
    BootstrappingKey(const SecretKey<Parameter>& s, Gen& rng) {
        for (std::size_t i = 0; i < Parameter::n; ++i) {
            data[i] = TRGSW<Parameter>::encrypt(s, (bool)s.get_lvl0(i), rng);
        }
    }
    template <RandGen Gen>
    static std::shared_ptr<BootstrappingKey> make_ptr(const SecretKey<Parameter>& s, Gen& rng) {
        return std::make_shared<BootstrappingKey>(s, rng);
    }
    const TRGSW<Parameter>& operator[](std::size_t i) const noexcept {
        return data[i];
    }
};

template <class Parameter>
void blind_rotate(TRLWE<Parameter>& trlwe, const TLWElvl0<Parameter>& tlwe, const BootstrappingKey<Parameter>& bk, const TRLWE<Parameter>& testvec) {
    static constexpr std::size_t bit_width = std::numeric_limits<Torus16>::digits;
    static constexpr std::size_t mask = (Parameter::N << 1) - 1;
    std::size_t b_tilda = (2 * Parameter::N - (tlwe.b() >> (bit_width - Parameter::Nbit - 1))) & mask;
    testvec.shift(trlwe, b_tilda);
    for (std::size_t i = 0; i < Parameter::n; ++i) {
        std::size_t a_tilda = ((tlwe.a(i) + static_cast<Torus16>(1 << (bit_width - 2 - Parameter::Nbit))) >> (bit_width - Parameter::Nbit - 1)) & mask;
        const TRLWE<Parameter> trlwe0 = trlwe;
        TRLWE<Parameter> trlwe1;
        trlwe0.shift(trlwe1, a_tilda);
        cmux(trlwe, bk[i], trlwe1, trlwe0);
    }
}

template <class Parameter>
constexpr TRLWE<Parameter> test_vector() {
    constexpr std::size_t bit_width = std::numeric_limits<Torus>::digits;
    constexpr Torus mu = 1 << (bit_width - 3);
    Poly<Torus, Parameter::N> m;
    for (auto& e : m) e = mu;
    return TRLWE<Parameter>::trivial_encrypt(m);
}

template <class Parameter>
void gate_bootstrapping_tlwe_to_tlwe(TLWElvl1<Parameter>& tlwe1, const TLWElvl0<Parameter>& tlwe, const BootstrappingKey<Parameter>& bk) {
    TRLWE<Parameter> trlwe;
    blind_rotate(trlwe, tlwe, bk, test_vector<Parameter>());
    sample_extract_index(tlwe1, trlwe, 0);
}