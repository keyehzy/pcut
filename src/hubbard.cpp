#include <pcut/models.hpp>
#include <bit>
#include <cmath>
#include <stdexcept>

namespace pcut::models {
LocalSpace hubbard_site() {
    return {{0,0,0,1},0,"Hubbard",{0,1,1,2},{0,1,1,0}};
}
Matrix hubbard_hopping(double amplitude) {
    if (!std::isfinite(amplitude)) throw std::invalid_argument("nonfinite Hubbard hopping");
    Matrix h=Matrix::Zero(16,16);
    for (unsigned state=0; state<16; ++state) for (unsigned spin=0; spin<2; ++spin)
        for (unsigned direction=0; direction<2; ++direction) {
            const unsigned from=2*direction+spin, to=2*(1-direction)+spin;
            if (!(state & (1u<<from)) || (state & (1u<<to))) continue;
            const unsigned removed=state^(1u<<from);
            const unsigned parity=static_cast<unsigned>(std::popcount(state & ((1u<<from)-1))+
                                                        std::popcount(removed & ((1u<<to)-1)))%2;
            h(removed|(1u<<to),state)+=amplitude*(parity ? -1 : 1);
        }
    return h;
}
PeriodicLattice hubbard_chain(double amplitude) {
    if (!std::isfinite(amplitude)) throw std::invalid_argument("nonfinite Hubbard hopping");
    return {1,{hubbard_site()},{{{{{0},0},{{1},0}},{{{hubbard_hopping(1),true},amplitude}}}},1};
}
PeriodicLattice hubbard_square(double amplitude) {
    if (!std::isfinite(amplitude)) throw std::invalid_argument("nonfinite Hubbard hopping");
    PeriodicLattice lattice{2,{hubbard_site()},{},1};
    for (const Coordinate& offset : {Coordinate{1,0},Coordinate{0,1}})
        lattice.interactions.push_back({{{{0,0},0},{offset,0}},{{{hubbard_hopping(1),true},amplitude}}});
    return lattice;
}
}
