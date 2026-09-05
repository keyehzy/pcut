#include <pcut/pcut.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Eigen/Eigenvalues>
#include <cmath>

TEST_CASE("Complex two-level model has the analytic pCUT eigenvalue series", "[model]") {
    pcut::Matrix v(2,2); v << 0,pcut::Complex(0,-2),pcut::Complex(0,2),0;
    const pcut::ClusterModel model({{{0,1},-0.7,"two-level"}},{{{0},v}},3.0);
    const pcut::EffectiveOperator effective(pcut::Coefficients({-1,0,1},6));
    const auto e = effective.vacuum(model);
    REQUIRE(e[0].real() == Catch::Approx(-0.7));
    REQUIRE(e[2].real() == Catch::Approx(-4.0/3));
    REQUIRE(e[4].real() == Catch::Approx(16.0/27));
    REQUIRE(e[6].real() == Catch::Approx(-128.0/243));
    const double x = 0.05;
    const double exact = -0.7+(3-std::sqrt(9+16*x*x))/2;
    REQUIRE(std::abs(pcut::evaluate(e,x).real()-exact) < 1e-10);
    const auto h = effective.block(model,{0,1});
    for (const auto& a : h) REQUIRE(a.isApprox(a.adjoint(),1e-12));
    for (unsigned n = 1; n <= 6; ++n) REQUIRE(std::abs(h[n](0,1)) < 1e-14);
    REQUIRE(h[2](1,1).real() == Catch::Approx(4.0/3));
}
TEST_CASE("Tensor roles, multi-site interactions and non-unit local charges", "[model]") {
    pcut::Matrix v = pcut::Matrix::Zero(12,12);
    v(11,0) = pcut::Complex(1,2); v(0,11) = std::conj(v(11,0));
    const pcut::ClusterModel model({{{0,1},0,"a"},{{0,2,3},0,"b"},{{0,1},0,"c"}},{{{2,1,0},v}});
    REQUIRE(model.encode({1,2,1}) == 11);
    REQUIRE(model.decode(11) == std::vector<unsigned>{1,2,1});
    REQUIRE(model.charge(11) == 5);
    REQUIRE(model.apply(5,{{0,1}}).at(11) == pcut::Complex(1,2));
    const pcut::EffectiveOperator effective(pcut::Coefficients({-5,5},2));
    REQUIRE(effective.vacuum(model)[2].real() == Catch::Approx(-1));
    REQUIRE_THROWS_AS(model.apply(5,{{0,1}},0),std::length_error);
    REQUIRE_THROWS_AS(model.encode({2,0,0}),std::out_of_range);
}
TEST_CASE("Finite dimer spectrum agrees with independent exact diagonalization", "[model][ed]") {
    const auto lattice = pcut::models::dimerized_chain(0.23);
    const auto model = pcut::cluster_model(lattice,{{0,{0}},{0,{1}}});
    const pcut::EffectiveOperator effective(pcut::Coefficients(pcut::charge_changes(lattice),4));
    const auto e = effective.vacuum(model);
    double previous = 0;
    for (double x : {0.12,0.06}) {
        const Eigen::SelfAdjointEigenSolver<pcut::Matrix> ed(model.dense_hamiltonian(x));
        REQUIRE(ed.info() == Eigen::Success);
        const double error = std::abs(pcut::evaluate(e,x).real()-ed.eigenvalues()[0]);
        REQUIRE(error < 1e-6);
        if (previous > 0) REQUIRE(previous/error > 24);
        previous = error;
    }
    REQUIRE_THROWS_AS(model.dense_hamiltonian(0.1,16),std::length_error);
}
TEST_CASE("Model contracts reject malformed input", "[model]") {
    pcut::Matrix x(2,2); x << 0,1,1,0;
    REQUIRE_THROWS_AS(pcut::ClusterModel({{{0,-1},0,"bad"}},{}),std::invalid_argument);
    REQUIRE_THROWS_AS(pcut::ClusterModel({{{0,1},0,"ok"}},{{{0,0},x}}),std::invalid_argument);
    REQUIRE_THROWS_AS(pcut::ClusterModel({{{0,1},0,"ok"}},{{{0},x}},0),std::invalid_argument);
    x(0,1)=3;
    REQUIRE_THROWS_AS(pcut::ClusterModel({{{0,1},0,"ok"}},{{{0},x}}),std::invalid_argument);
    const pcut::ClusterModel model({{{0,1},0,"ok"}},{});
    REQUIRE_THROWS_AS(pcut::EffectiveOperator(pcut::Coefficients({0},1)).block(model,{0,0}),std::invalid_argument);
    REQUIRE_THROWS_AS(pcut::ClusterModel(std::vector<pcut::LocalSpace>(64,{{0,1},0,"spin"}),{}),std::length_error);
}

TEST_CASE("Physical energy unit rescaling preserves the effective series", "[model][units]") {
    const pcut::EffectiveOperator effective(pcut::Coefficients({-1,1},6));
    for (double scale : {1e-150,1e-90,1.0,1e90,1e150}) {
        CAPTURE(scale);
        pcut::Matrix x(2,2); x<<0,pcut::Complex(0,-scale),pcut::Complex(0,scale),0;
        const pcut::ClusterModel model({{{0,1},0,"two-level"}},{{{0},x}},scale);
        const auto block=effective.block(model,{0,1});
        const auto vacuum=effective.vacuum(model);
        for (unsigned n : {2U,4U,6U}) {
            const double expected=n==2 ? -1 : n==4 ? 1 : -2;
            REQUIRE((vacuum[n]/scale).real()==Catch::Approx(expected));
            REQUIRE((block[n](0,0)/scale).real()==Catch::Approx(expected));
            REQUIRE((block[n](1,1)/scale).real()==Catch::Approx(-expected));
        }
    }
}
TEST_CASE("Dense block budgets count every perturbative order", "[model][budget]") {
    const pcut::ClusterModel model({{{0,1},0,"spin"}},{});
    const pcut::EffectiveOperator effective(pcut::Coefficients({0},64));
    pcut::SolverOptions options; options.max_matrix_elements=259;
    REQUIRE_THROWS_AS(effective.block(model,{0,1},options),std::length_error);
    options.max_matrix_elements=260;
    REQUIRE(effective.block(model,{0,1},options).size()==65);
}
