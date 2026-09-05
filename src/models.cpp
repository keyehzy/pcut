#include <pcut/models.hpp>
#include <unsupported/Eigen/KroneckerProduct>
#include <array>
#include <cmath>
#include <stdexcept>

namespace pcut::models {
namespace {
// Tensor index is first factor + dim(first)*second factor.
Matrix tensor(const Matrix& first, const Matrix& second) {
    return Eigen::kroneckerProduct(second,first).eval();
}
std::array<std::array<Matrix,3>,2> dimer_spins() {
    const Complex i(0,1);
    Matrix sx(2,2), sy(2,2), sz(2,2);
    sx << 0,0.5,0.5,0;
    sy << 0,-0.5*i,0.5*i,0;
    sz << 0.5,0,0,-0.5;
    const Matrix identity = Matrix::Identity(2,2);
    Matrix u = Matrix::Zero(4,4);
    const double r = 1/std::sqrt(2.0);
    u(2,0)=r; u(1,0)=-r; // |up,down> - |down,up>
    u(0,1)=1; u(2,2)=r; u(1,2)=r; u(3,3)=1;
    std::array<std::array<Matrix,3>,2> spins;
    const std::array<Matrix,3> axes{sx,sy,sz};
    for (std::size_t a = 0; a < 3; ++a) {
        spins[0][a] = u.adjoint()*tensor(axes[a],identity)*u;
        spins[1][a] = u.adjoint()*tensor(identity,axes[a])*u;
    }
    return spins;
}
}
LocalSpace spin_dimer(double J) {
    if (!std::isfinite(J) || J <= 0) throw std::invalid_argument("J must be positive and finite");
    return {{0,1,1,1},-0.75*J,"spin-1/2 dimer",{},{} };
}
Matrix dimer_bond(double alpha, double J) {
    if (!std::isfinite(alpha)) throw std::invalid_argument("alpha must be finite");
    (void)spin_dimer(J);
    const auto spin = dimer_spins();
    Matrix v = Matrix::Zero(16,16);
    for (std::size_t a = 0; a < 3; ++a)
        v += J*(tensor(spin[1][a],spin[0][a]) + alpha*(tensor(spin[0][a],spin[0][a])+tensor(spin[1][a],spin[1][a])));
    // Clean only roundoff from this analytic basis transformation, not user matrices.
    for (Eigen::Index i = 0; i < v.size(); ++i) if (std::abs(v.data()[i]) < 1e-14*J) v.data()[i] = 0;
    return v;
}
PeriodicLattice dimerized_chain(double alpha, double J) {
    if (!std::isfinite(alpha)) throw std::invalid_argument("alpha must be finite");
    const Matrix base=dimer_bond(0,J);
    const Matrix frustration=dimer_bond(1,J)-base;
    return {1,{spin_dimer(J)},{{{{{0},0},{{1},0}},{{{base,false},1},{{frustration,false},alpha}}}},J};
}
PeriodicLattice ising_chain(double coupling) {
    if (!std::isfinite(coupling)) throw std::invalid_argument("coupling must be finite");
    Matrix x(2,2); x << 0,1,1,0;
    return {1,{{{0,1},0,"hardcore spin",{},{} }},{{{{{0},0},{{1},0}},{{{-tensor(x,x),false},coupling}}}},1.0};
}
}
