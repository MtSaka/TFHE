#pragma once

#include "tlwe.cpp"
#include "bootstrapping.cpp"

template <class Parameter>
TLWElvl0<Parameter> hom_nand(const TLWElvl0<Parameter>& x, const TLWElvl0<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus16 mu = (1 << (std::numeric_limits<Torus16>::digits - 3));
    TLWElvl0<Parameter> tlwe0 = {};
    tlwe0.b() = mu;
    tlwe0 -= x + y;
    TLWElvl1<Parameter> tlwe1 = gate_bootstrapping_tlwe_to_tlwe(tlwe0, bk);
    return identity_key_switch(tlwe1, ks);
}