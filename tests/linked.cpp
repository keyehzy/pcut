#include <pcut/pcut.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

TEST_CASE("Chain catalog includes every embedded subinterval", "[lattice]") {
    const pcut::ClusterCatalog catalog(pcut::models::ising_chain(),5);
    REQUIRE(catalog.entries().size() == 5);
    for (std::size_t n = 1; n <= 5; ++n) {
        const auto& e = catalog.entries()[n-1];
        REQUIRE(e.sites.size() == n+1);
        REQUIRE(e.subclusters.size() == n*(n+1)/2-1);
        REQUIRE(pcut::connected(catalog.lattice(),e.edges));
    }
    const auto a = pcut::normalize({{0,{-7}},{0,{-6}}});
    REQUIRE(a.cluster == catalog.entries()[1].edges);
    REQUIRE(a.shift == pcut::Coordinate{-7});
}
TEST_CASE("Square lattice counts colored orientations and loops", "[lattice]") {
    auto lattice = pcut::models::ising_chain();
    const auto v = lattice.interactions[0].matrix;
    lattice.dimension=2;
    lattice.interactions={{{{{0,0},0},{{1,0},0}},v},{{{{0,0},0},{{0,1},0}},v}};
    const pcut::ClusterCatalog catalog(lattice,4);
    std::vector<unsigned> counts(5);
    unsigned loops = 0;
    for (const auto& e : catalog.entries()) {
        ++counts[e.edges.size()];
        if (e.edges.size()==4 && e.sites.size()==4) ++loops;
    }
    REQUIRE(counts[1]==2);
    REQUIRE(counts[2]==6);
    REQUIRE(counts[3]==22);
    REQUIRE(counts[4]==88);
    REQUIRE(loops==1);
    REQUIRE_THROWS_AS(pcut::ClusterCatalog(lattice,4,{3,100}),std::length_error);
}
TEST_CASE("Scalar subtraction cancels an extensive bond contribution", "[linked]") {
    const pcut::ClusterCatalog catalog(pcut::models::ising_chain(),4);
    const auto result = pcut::linked_scalar(catalog,4,[](const pcut::ClusterModel& m) {
        pcut::Series s(5); s[1]=2.0*static_cast<double>(m.sites()-1); return s;
    },-0.5);
    REQUIRE(result.per_cell[0].real()==Catch::Approx(-0.5));
    REQUIRE(result.per_cell[1].real()==Catch::Approx(2));
    for (std::size_t i=1;i<result.weights.size();++i) REQUIRE(std::abs(result.weights[i][1])<1e-14);
}
TEST_CASE("Ising thermodynamic energy and dispersion match the free-fermion expansion", "[linked][analytic]") {
    const auto lattice=pcut::models::ising_chain();
    const pcut::ClusterCatalog catalog(lattice,4);
    const pcut::EffectiveOperator effective(pcut::Coefficients(pcut::charge_changes(lattice),4));
    const auto result=pcut::linked_expand(catalog,effective);
    REQUIRE(result.energy_per_cell[2].real()==Catch::Approx(-0.5));
    REQUIRE(result.energy_per_cell[4].real()==Catch::Approx(-0.125));
    for (double k : {0.0,0.4,1.7,3.141592653589793}) {
        const auto w=result.bloch_series({k});
        const double c=std::cos(k);
        REQUIRE(w[0](0,0).real()==Catch::Approx(1));
        REQUIRE(w[1](0,0).real()==Catch::Approx(-2*c));
        REQUIRE(w[2](0,0).real()==Catch::Approx(2-2*c*c).margin(1e-12));
        REQUIRE(w[3](0,0).real()==Catch::Approx(4*c-4*c*c*c).margin(1e-12));
        REQUIRE(w[4](0,0).real()==Catch::Approx(-2+12*c*c-10*std::pow(c,4)).margin(1e-12));
    }
    REQUIRE_THROWS_AS(pcut::linked_expand(pcut::ClusterCatalog(lattice,2),effective),std::invalid_argument);
    REQUIRE_THROWS_AS(result.bloch({0,0},0.1),std::invalid_argument);
}
TEST_CASE("On-site terms, multiple basis sites and complex directed hopping", "[linked]") {
    pcut::Matrix onsite(2,2); onsite<<0,0,0,2;
    pcut::Matrix bond=pcut::Matrix::Zero(4,4);
    bond(2,1)=pcut::Complex(0,1); bond(1,2)=pcut::Complex(0,-1);
    pcut::PeriodicLattice lattice{1,{{{0,1},-0.2,"a"},{{0,1},-0.3,"b"}},
        {{{{{0},0}},onsite},{{{{0},0},{{0},1}},bond},{{{{0},1},{{1},0}},bond}},1};
    const pcut::EffectiveOperator effective(pcut::Coefficients(pcut::charge_changes(lattice),3));
    const auto result=pcut::linked_expand(pcut::ClusterCatalog(lattice,3),effective);
    REQUIRE(result.energy_per_cell[0].real()==Catch::Approx(-0.5));
    const auto h=result.bloch_series({0.7});
    REQUIRE(h[1].isApprox(h[1].adjoint(),1e-12));
    REQUIRE(h[1](0,0).real()==Catch::Approx(2));
    const pcut::Complex expected=pcut::Complex(0,-1)+pcut::Complex(0,1)*std::polar(1.0,-0.7);
    REQUIRE(std::abs(h[1](0,1)-expected)<1e-12);
    REQUIRE(h[2].norm()<1e-12);
    REQUIRE(h[3].norm()<1e-12);
}
