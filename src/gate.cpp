#pragma once

#include "tlwe.cpp"
#include "bootstrapping.cpp"

template <class Parameter>
void hom_nand(TLWElvl0<Parameter>& res, const TLWElvl0<Parameter>& x, const TLWElvl0<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus16 mu = (1 << (std::numeric_limits<Torus16>::digits - 3));
    TLWElvl0<Parameter> tlwe0 = -x - y + mu;
    TLWElvl1<Parameter> tlwe1;
    gate_bootstrapping_tlwe_to_tlwe(tlwe1, tlwe0, bk);
    identity_key_switch(res, tlwe1, ks);
}

template <class Parameter>
void hom_and(TLWElvl0<Parameter>& res, const TLWElvl0<Parameter>& x, const TLWElvl0<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus16 mu = (1 << (std::numeric_limits<Torus16>::digits - 3));
    TLWElvl0<Parameter> tlwe0 = x + y - mu;
    TLWElvl1<Parameter> tlwe1;
    gate_bootstrapping_tlwe_to_tlwe(tlwe1, tlwe0, bk);
    identity_key_switch(res, tlwe1, ks);
}

template <class Parameter>
void hom_or(TLWElvl0<Parameter>& res, const TLWElvl0<Parameter>& x, const TLWElvl0<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus16 mu = (1 << (std::numeric_limits<Torus16>::digits - 3));
    TLWElvl0<Parameter> tlwe0 = x + y + mu;
    TLWElvl1<Parameter> tlwe1;
    gate_bootstrapping_tlwe_to_tlwe(tlwe1, tlwe0, bk);
    identity_key_switch(res, tlwe1, ks);
}

template <class Parameter>
void hom_xor(TLWElvl0<Parameter>& res, const TLWElvl0<Parameter>& x, const TLWElvl0<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus16 mu = (1 << (std::numeric_limits<Torus16>::digits - 3));
    TLWElvl0<Parameter> tlwe0 = (x + y + mu);
    tlwe0 *= 2;
    TLWElvl1<Parameter> tlwe1;
    gate_bootstrapping_tlwe_to_tlwe(tlwe1, tlwe0, bk);
    identity_key_switch(res, tlwe1, ks);
}
template <class Parameter>
void hom_nor(TLWElvl0<Parameter>& res, const TLWElvl0<Parameter>& x, const TLWElvl0<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus16 mu = (1 << (std::numeric_limits<Torus16>::digits - 3));
    TLWElvl0<Parameter> tlwe0 = -x - y - mu;
    TLWElvl1<Parameter> tlwe1;
    gate_bootstrapping_tlwe_to_tlwe(tlwe1, tlwe0, bk);
    identity_key_switch(res, tlwe1, ks);
}
template <class Parameter>
void hom_xnor(TLWElvl0<Parameter>& res, const TLWElvl0<Parameter>& x, const TLWElvl0<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus16 mu = (1 << (std::numeric_limits<Torus16>::digits - 3));
    TLWElvl0<Parameter> tlwe0 = (x - y + mu);
    tlwe0 *= 2;
    TLWElvl1<Parameter> tlwe1;
    gate_bootstrapping_tlwe_to_tlwe(tlwe1, tlwe0, bk);
    identity_key_switch(res, tlwe1, ks);
}

template <class Parameter>
void hom_mux(TLWElvl0<Parameter>& res, const TLWElvl0<Parameter>& x, const TLWElvl0<Parameter>& y, const TLWElvl0<Parameter>& s, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    static constexpr Torus16 mu = (1 << (std::numeric_limits<Torus16>::digits - 3));
    TLWElvl0<Parameter> xs, ys;
    hom_and(xs, x, s, bk, ks);
    hom_and(ys, y, -s, bk, ks);
    res = xs + ys + mu;
    // hom_or(res, xs, ys, bk, ks);
    /*
    TLWElvl0 tlwe0 = xs + ys + mu;
    TLWElvl1<Parameter> tlwe1;
    gate_bootstrapping_tlwe_to_tlwe(tlwe1, tlwe0, bk);
    identity_key_switch(res, tlwe1, ks);*/
}

template <class Parameter>
void hom_half_adder(TLWElvl0<Parameter>& ress, TLWElvl0<Parameter>& resc, const TLWElvl0<Parameter>& x, const TLWElvl0<Parameter>& y, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    hom_xor(ress, x, y, bk, ks);
    hom_and(resc, x, y, bk, ks);
}
template <class Parameter>
void hom_full_adder(TLWElvl0<Parameter>& ress, TLWElvl0<Parameter>& resc, const TLWElvl0<Parameter>& x, const TLWElvl0<Parameter>& y, const TLWElvl0<Parameter>& c, const BootstrappingKey<Parameter>& bk, const KeySwitchKey<Parameter>& ks) {
    TLWElvl0<Parameter> xxory;
    hom_xor(xxory, x, y, bk, ks);
    TLWElvl0<Parameter> s;
    hom_xor(ress, xxory, c, bk, ks);
    TLWElvl0<Parameter> t1, t2;
    hom_and(t1, xxory, c, bk, ks);
    hom_and(t2, x, y, bk, ks);
    hom_or(resc, t1, t2, bk, ks);
}