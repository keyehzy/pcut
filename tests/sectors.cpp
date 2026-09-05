#include <pcut/pcut.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

TEST_CASE("General sector linking agrees with the specialized vacuum and hopping driver", "[sectors]") {
    const auto lattice=pcut::models::dimerized_chain(0.17);
    const pcut::ClusterCatalog catalog(lattice,4);
    const pcut::EffectiveOperator effective(pcut::Coefficients(pcut::charge_changes(lattice),4));
    const auto specialized=pcut::linked_expand(catalog,effective);
    const auto general=pcut::linked_expand_sectors(catalog,effective,1);
    for (unsigned n=0;n<=4;++n) REQUIRE(std::abs(specialized.energy_per_cell[n]-general.energy_per_cell[n])<1e-12);
    for (const auto& [key,series] : general.kernels) {
        REQUIRE(key.output.size()==1);
        REQUIRE(key.input.size()==1);
        const pcut::Hopping h{key.output[0].local-1,key.input[0].local-1,
            {key.output[0].site.cell[0]-key.input[0].site.cell[0]}};
        const auto& expected=specialized.hopping.at(h);
        for (unsigned n=0;n<=4;++n) REQUIRE(std::abs(series[n]-expected[n])<1e-12);
    }
}
TEST_CASE("Three-body kernels retain the interaction after spectator subtraction", "[sectors][hypergraph]") {
    pcut::Matrix v=pcut::Matrix::Zero(8,8); v(7,7)=2.3;
    pcut::PeriodicLattice lattice{1,{{{0,1},-0.2,"spin"}},{{{{{0},0},{{1},0},{{2},0}},v}},1};
    const pcut::ClusterCatalog catalog(lattice,3);
    const pcut::EffectiveOperator effective(pcut::Coefficients({0},3));
    const auto result=pcut::linked_expand_sectors(catalog,effective,3);
    unsigned triples=0;
    for (const auto& [key,s] : result.kernels) {
        if (key.output.size()==3) {
            REQUIRE(key.output==key.input);
            REQUIRE(s[1].real()==Catch::Approx(2.3)); ++triples;
        } else REQUIRE(std::abs(s[1])<1e-12);
        REQUIRE(std::abs(s[2])<1e-12);
        REQUIRE(std::abs(s[3])<1e-12);
    }
    REQUIRE(triples==1);
    REQUIRE(result.energy_per_cell[0].real()==Catch::Approx(-0.2));
}
TEST_CASE("Irreducible kernels remove disconnected dressed spectators", "[sectors]") {
    pcut::Matrix x(2,2); x<<0,1,1,0;
    const pcut::ClusterModel model({{{0,1},-0.4,"a"},{{0,1},-0.4,"b"}},{{{0},x},{{1},x}});
    const pcut::EffectiveOperator effective(pcut::Coefficients({-1,0,1},4));
    const auto result=pcut::irreducible_sectors(model,effective,2);
    REQUIRE(result.basis.back()==3);
    for (const auto& k : result.kernels) REQUIRE(std::abs(k(3,3))<1e-12);
    REQUIRE(result.kernels[2](1,1).real()==Catch::Approx(2));
    REQUIRE(result.kernels[4](1,1).real()==Catch::Approx(-2));
    REQUIRE_THROWS_AS(pcut::irreducible_sectors(model,effective,2,{2,100,100,{}}),std::length_error);
}
TEST_CASE("Unequal local charges permit a one-to-two excitation conversion kernel", "[sectors]") {
    pcut::Matrix v=pcut::Matrix::Zero(8,8);
    v(1,6)=0.7; v(6,1)=0.7;
    const pcut::ClusterModel model({{{0,2},0,"q2"},{{0,1},0,"q1"},{{0,1},0,"q1"}},{{{0,1,2},v}});
    const pcut::EffectiveOperator effective(pcut::Coefficients({0},2));
    const auto result=pcut::irreducible_sectors(model,effective,2);
    const auto i=static_cast<Eigen::Index>(std::find(result.basis.begin(),result.basis.end(),1)-result.basis.begin());
    const auto j=static_cast<Eigen::Index>(std::find(result.basis.begin(),result.basis.end(),6)-result.basis.begin());
    REQUIRE(result.kernels[1](i,j).real()==Catch::Approx(0.7));
    REQUIRE(result.kernels[2].norm()<1e-12);
}

TEST_CASE("Ising two-particle kernels recover the Jordan-Wigner correlated hopping", "[sectors][analytic]") {
    const auto lattice=pcut::models::ising_chain();
    const pcut::EffectiveOperator effective(pcut::Coefficients(pcut::charge_changes(lattice),4));
    const auto result=pcut::linked_expand_sectors(pcut::ClusterCatalog(lattice,4),effective,2);
    const pcut::Excitation left{{{0},0},1}, middle{{{1},0},1}, right{{{2},0},1};
    const auto& correlated=result.kernels.at({{middle,right},{left,middle}});
    // Free-fermion t_2=-lambda^2/2 times the string (1-2 n_middle).
    REQUIRE(correlated[2].real()==Catch::Approx(1));
    const auto& hopping=result.kernels.at({{right},{left}});
    REQUIRE(hopping[2].real()==Catch::Approx(-0.5));
}

TEST_CASE("Sector storage is cumulative and excess catalog orders are skipped", "[sectors][budget]") {
    const auto lattice=pcut::models::ising_chain();
    const pcut::EffectiveOperator first(pcut::Coefficients(pcut::charge_changes(lattice),1));
    pcut::SectorOptions options; options.max_matrix_elements=30;
    const auto expected=pcut::linked_expand_sectors(pcut::ClusterCatalog(lattice,1),first,1,options);
    const auto actual=pcut::linked_expand_sectors(pcut::ClusterCatalog(lattice,3),first,1,options);
    REQUIRE(actual.energy_per_cell==expected.energy_per_cell);
    REQUIRE(actual.kernels==expected.kernels);
    const pcut::EffectiveOperator second(pcut::Coefficients(pcut::charge_changes(lattice),2));
    options.solver.max_matrix_elements=48;
    options.max_matrix_elements=70; // both individual blocks fit, their sum does not
    REQUIRE_THROWS_AS(pcut::linked_expand_sectors(pcut::ClusterCatalog(lattice,2),second,1,options),std::length_error);
    options.max_matrix_elements=26;
    REQUIRE_THROWS_AS(pcut::irreducible_sectors(pcut::cluster_model(lattice,{{0,{0}}}),second,1,options),std::length_error);
}
