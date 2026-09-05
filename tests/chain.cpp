#include <pcut/pcut.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

TEST_CASE("Dimer chain vacuum matches literature through sixth order", "[chain][literature]") {
    const pcut::EffectiveOperator effective(pcut::Coefficients({-2,-1,0,1,2},6));
    for (double a : {0.0,0.17,0.5}) {
        INFO("alpha=" << a);
        const auto result=pcut::linked_expand(pcut::ClusterCatalog(pcut::models::dimerized_chain(a),6),effective,{false});
        const double b=std::pow(1-2*a,2);
        // cond-mat/9906243 Eq. E_grund: per SPIN, lambda_bar=lambda/4;
        // restore the physical -3/8 reference and multiply by two per dimer.
        const std::vector<double> expected{-0.75,0,
            -2*b*(3.0/4)/16,
            -2*b*(3.0/4+3.0*a/2)/64,
            -2*b*(13.0/16+27.0*a/4-3.0*a*a/4)/256,
            -2*b*(89.0/48+311.0*a/24+93.0*a*a/4-45.0*a*a*a/2)/1024,
            -2*b*(463.0/96+227.0*a/9+1307.0*a*a/12-42*a*a*a-159.0*std::pow(a,4)/2)/4096};
        for (unsigned n=0;n<=6;++n) {
            INFO("order=" << n);
            REQUIRE(result.energy_per_cell[n].real()==Catch::Approx(expected[n]).margin(2e-11));
            REQUIRE(std::abs(result.energy_per_cell[n].imag())<1e-13);
        }
    }
}
TEST_CASE("Dimer one-particle hopping obeys first-order spin algebra", "[chain]") {
    for (double a : {0.0,0.23,0.5}) {
        const auto lattice=pcut::models::dimerized_chain(a);
        const pcut::EffectiveOperator effective(pcut::Coefficients({-2,-1,0,1,2},3));
        const auto result=pcut::linked_expand(pcut::ClusterCatalog(lattice,3),effective);
        const auto it=result.hopping.find({0,0,{1}});
        REQUIRE(it!=result.hopping.end());
        REQUIRE(it->second[1].real()==Catch::Approx(-(1-2*a)/4).margin(1e-12));
        const auto w=result.bloch_series({0.41});
        for (const auto& h : w) {
            REQUIRE(h.isApprox(h.adjoint(),1e-12));
            REQUIRE(std::abs(h(0,0)-h(1,1))<1e-12);
            REQUIRE(std::abs(h(0,0)-h(2,2))<1e-12);
        }
    }
}

TEST_CASE("Dimer dispersion matches Appendix D through sixth order", "[chain][literature]") {
    const pcut::EffectiveOperator effective(pcut::Coefficients({-2,-1,0,1,2},6));
    for (double a : {0.0,0.23,0.5}) {
        INFO("alpha=" << a);
        const auto result=pcut::linked_expand(pcut::ClusterCatalog(pcut::models::dimerized_chain(a),6),effective);
        std::ifstream file(PCUT_TEST_DATA "/chain_hopping.tsv");
        REQUIRE(file.good());
        std::string line;
        unsigned checked=0;
        while (std::getline(file,line)) {
            if (line.empty() || line[0]=='#') continue;
            std::istringstream input(line);
            int r; unsigned n; std::string polynomial;
            input>>r>>n>>polynomial;
            std::replace(polynomial.begin(),polynomial.end(),',',' ');
            std::istringstream numbers(polynomial);
            std::string number;
            double expected=0, power=1;
            while (numbers>>number) {
                expected+=pcut::Rational(number).convert_to<double>()*power;
                power*=1-2*a;
            }
            expected/=std::pow(4.0,n);
            if (r>0) expected/=2; // appendix gives cosine amplitudes, verified independently at order 1
            INFO("r=" << r << " order=" << n);
            for (std::size_t flavor=0;flavor<3;++flavor) {
                const auto it=result.hopping.find({flavor,flavor,{r}});
                REQUIRE(it!=result.hopping.end());
                REQUIRE(it->second[n].real()==Catch::Approx(expected).margin(5e-11));
            }
            ++checked;
        }
        REQUIRE(checked==27);
    }
}
