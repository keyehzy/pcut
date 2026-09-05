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

TEST_CASE("Linked scalar and particle drivers skip catalog orders they do not need", "[linked]") {
    const auto lattice=pcut::models::ising_chain();
    const pcut::ClusterCatalog small(lattice,1), large(lattice,3);
    const pcut::EffectiveOperator effective(pcut::Coefficients(pcut::charge_changes(lattice),1));
    const auto expected=pcut::linked_expand(small,effective);
    const auto actual=pcut::linked_expand(large,effective);
    REQUIRE(actual.energy_per_cell==expected.energy_per_cell);
    REQUIRE(actual.hopping==expected.hopping);
    REQUIRE(actual.vacuum_weights.size()==1);
    unsigned calls=0;
    const auto scalar=pcut::linked_scalar(large,1,[&](const auto& model) {
        ++calls; REQUIRE(model.sites()==2); return pcut::Series{0,2};
    });
    REQUIRE(calls==1);
    REQUIRE(scalar.weights.size()==1);
    REQUIRE(scalar.per_cell[1]==pcut::Complex(2));
    const auto zero=pcut::linked_scalar(large,0,[&](const auto&) { ++calls; return pcut::Series{0}; });
    REQUIRE(calls==1);
    REQUIRE(zero.weights.empty());
}
TEST_CASE("Coupling sweeps bind shared topology and preserve ordered compiled legs", "[lattice][sweep]") {
    auto lattice=pcut::models::dimerized_chain(0.17);
    const pcut::ClusterCatalog original(lattice,2);
    lattice=pcut::models::dimerized_chain(0.31);
    lattice.gap=2.5;
    lattice.cell[0].vacuum_energy=-0.9;
    const pcut::ClusterCatalog rebound(lattice,original.topology()), fresh(lattice,2);
    REQUIRE(&rebound.entries()==&original.entries());
    const pcut::EffectiveOperator effective(pcut::Coefficients(pcut::charge_changes(lattice),2));
    const auto a=pcut::linked_expand(rebound,effective), b=pcut::linked_expand(fresh,effective);
    REQUIRE(a.energy_per_cell==b.energy_per_cell);
    REQUIRE(a.hopping==b.hopping);
    for (const auto& entry : rebound.entries()) {
        const auto compiled=rebound.model(entry.edges).dense_hamiltonian(0.3);
        REQUIRE(compiled.isApprox(pcut::cluster_model(lattice,entry.edges).dense_hamiltonian(0.3),1e-13));
    }
    const auto surviving_model=[&] {
        const pcut::ClusterCatalog temporary(lattice,original.topology());
        return temporary.model({{0,{0}},{0,{1}}});
    }();
    REQUIRE(surviving_model.dense_hamiltonian(0.3).isApprox(
        pcut::cluster_model(lattice,{{0,{0}},{0,{1}}}).dense_hamiltonian(0.3),1e-13));
    auto zero=lattice; zero.interactions[0].matrix.setZero();
    const auto uncoupled=pcut::linked_expand(pcut::ClusterCatalog(zero,original.topology()),effective);
    REQUIRE(uncoupled.energy_per_cell[0]==pcut::Complex(-0.9));
    REQUIRE(uncoupled.energy_per_cell[2]==pcut::Complex{});
    auto bad=lattice; bad.interactions[0].legs[0].cell[0]+=2;
    REQUIRE_THROWS_AS(pcut::ClusterCatalog(bad,original.topology()),std::invalid_argument);
    bad=lattice; bad.interactions[0].matrix(0,1)+=1;
    REQUIRE_THROWS_AS(pcut::ClusterCatalog(bad,original.topology()),std::invalid_argument);

    pcut::Matrix v=pcut::Matrix::Zero(6,6);
    v(1,2)=pcut::Complex(0,0.7); v(2,1)=std::conj(v(1,2));
    pcut::PeriodicLattice mixed{1,{{{0,1},0,"a"},{{0,1,2},0,"b"}},
        {{{{{1},1},{{0},0}},v}},1};
    const pcut::ClusterCatalog colored(mixed,1);
    const pcut::Cluster edges{{0,{-2}},{0,{1}}};
    REQUIRE(colored.model(edges).dense_hamiltonian(0.2).isApprox(
        pcut::cluster_model(mixed,edges).dense_hamiltonian(0.2),1e-13));
}
