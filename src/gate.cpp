#pragma once

#include "tlwe.cpp"
#include "bootstrapping.cpp"

template <class Parameter>
void hom_nand(TLWElvl1<Parameter>& res, const TLWElvl1<Parameter>& x, const TLWElvl1<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus mu = (1 << (std::numeric_limits<Torus>::digits - 3));
    TLWElvl1<Parameter> tlwe1 = -x - y + mu;
    TLWElvl0<Parameter> tlwe0;
    identity_key_switch(tlwe0, tlwe1, ks);
    gate_bootstrapping_tlwe_to_tlwe(res, tlwe0, bk);
}

template <class Parameter>
void hom_and(TLWElvl1<Parameter>& res, const TLWElvl1<Parameter>& x, const TLWElvl1<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus mu = (1 << (std::numeric_limits<Torus>::digits - 3));
    TLWElvl1<Parameter> tlwe1 = x + y - mu;
    TLWElvl0<Parameter> tlwe0;
    identity_key_switch(tlwe0, tlwe1, ks);
    gate_bootstrapping_tlwe_to_tlwe(res, tlwe0, bk);
}

template <class Parameter>
void hom_or(TLWElvl1<Parameter>& res, const TLWElvl1<Parameter>& x, const TLWElvl1<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus mu = (1 << (std::numeric_limits<Torus>::digits - 3));
    TLWElvl1<Parameter> tlwe1 = x + y + mu;
    TLWElvl0<Parameter> tlwe0;
    identity_key_switch(tlwe0, tlwe1, ks);
    gate_bootstrapping_tlwe_to_tlwe(res, tlwe0, bk);
}

template <class Parameter>
void hom_xor(TLWElvl1<Parameter>& res, const TLWElvl1<Parameter>& x, const TLWElvl1<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus mu = (1 << (std::numeric_limits<Torus>::digits - 3));
    TLWElvl1<Parameter> tlwe1 = (x + y + mu);
    tlwe1 *= 2;
    TLWElvl0<Parameter> tlwe0;
    identity_key_switch(tlwe0, tlwe1, ks);
    gate_bootstrapping_tlwe_to_tlwe(res, tlwe0, bk);
}
template <class Parameter>
void hom_nor(TLWElvl1<Parameter>& res, const TLWElvl1<Parameter>& x, const TLWElvl1<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus mu = (1 << (std::numeric_limits<Torus>::digits - 3));
    TLWElvl1<Parameter> tlwe1 = -x - y - mu;
    TLWElvl0<Parameter> tlwe0;
    identity_key_switch(tlwe0, tlwe1, ks);
    gate_bootstrapping_tlwe_to_tlwe(res, tlwe0, bk);
}
template <class Parameter>
void hom_xnor(TLWElvl1<Parameter>& res, const TLWElvl1<Parameter>& x, const TLWElvl1<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus mu = (1 << (std::numeric_limits<Torus>::digits - 3));
    TLWElvl1<Parameter> tlwe1 = (x - y + mu);
    tlwe1 *= 2;
    TLWElvl0<Parameter> tlwe0;
    identity_key_switch(tlwe0, tlwe1, ks);
    gate_bootstrapping_tlwe_to_tlwe(res, tlwe0, bk);
}

template <class Parameter>
void hom_mux(TLWElvl1<Parameter>& res, const TLWElvl1<Parameter>& x, const TLWElvl1<Parameter>& y, const TLWElvl1<Parameter>& s, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus mu = (1 << (std::numeric_limits<Torus>::digits - 3));
    TLWElvl1<Parameter> xs, ys;
    hom_and(xs, x, s, bk, ks);
    hom_and(ys, y, -s, bk, ks);
    res = xs + ys + mu;
}

template <class Parameter>
void hom_half_adder(TLWElvl1<Parameter>& ress, TLWElvl1<Parameter>& resc, const TLWElvl1<Parameter>& x, const TLWElvl1<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    hom_xor(ress, x, y, bk, ks);
    hom_and(resc, x, y, bk, ks);
}
template <class Parameter>
void hom_full_adder(TLWElvl1<Parameter>& ress, TLWElvl1<Parameter>& resc, const TLWElvl1<Parameter>& x, const TLWElvl1<Parameter>& y, const TLWElvl1<Parameter>& c, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    TLWElvl1<Parameter> xxory;
    hom_xor(xxory, x, y, bk, ks);
    TLWElvl1<Parameter> s;
    hom_xor(ress, xxory, c, bk, ks);
    TLWElvl1<Parameter> t1, t2;
    hom_and(t1, xxory, c, bk, ks);
    hom_and(t2, x, y, bk, ks);
    hom_or(resc, t1, t2, bk, ks);
}