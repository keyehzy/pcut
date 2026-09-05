#include <pcut/pcut.hpp>
#include <Eigen/Eigenvalues>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        unsigned order = 4;
        double alpha = 0, lambda = 0.3, momentum = 0;
        bool particles = true;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "pcut_chain [--order N] [--alpha A] [--lambda L] [--k K] [--vacuum-only]\n"
                          << "J=1; energy per dimer includes -3/4. k uses dimer spacing.\n";
                return 0;
            }
            if (arg == "--vacuum-only") { particles = false; continue; }
            if (i+1 == argc) throw std::invalid_argument("missing value for " + arg);
            const std::string value = argv[++i];
            std::size_t end = 0;
            if (arg == "--order") {
                const auto parsed = std::stoll(value,&end);
                if (parsed < 0 || parsed > 64) throw std::invalid_argument("order must be between 0 and 64");
                order = static_cast<unsigned>(parsed);
            } else {
                const double parsed = std::stod(value,&end);
                if (!std::isfinite(parsed)) throw std::invalid_argument("arguments must be finite");
                if (arg == "--alpha") alpha = parsed;
                else if (arg == "--lambda") lambda = parsed;
                else if (arg == "--k") momentum = parsed;
                else throw std::invalid_argument("unknown option " + arg);
            }
            if (end != value.size()) throw std::invalid_argument("invalid number " + value);
        }
        const auto start = std::chrono::steady_clock::now();
        const auto lattice = pcut::models::dimerized_chain(alpha);
        const pcut::Coefficients coefficients(pcut::charge_changes(lattice),order);
        const pcut::EffectiveOperator effective(coefficients);
        const pcut::ClusterCatalog catalog(lattice,order);
        const auto result = pcut::linked_expand(catalog,effective,{particles,{}});
        std::cout << std::setprecision(15);
        std::cout << "# alpha=" << alpha << " order=" << order << " clusters=" << catalog.entries().size()
                  << " universal_terms=" << coefficients.terms().size() << '\n';
        std::cout << "order,energy_per_dimer,energy_per_spin";
        if (particles) std::cout << ",omega_k";
        std::cout << '\n';
        const auto dispersion = result.bloch_series({momentum});
        for (unsigned n = 0; n <= order; ++n) {
            std::cout << n << ',' << result.energy_per_cell[n].real() << ',' << result.energy_per_cell[n].real()/2;
            if (particles) std::cout << ',' << dispersion[n](0,0).real();
            std::cout << '\n';
        }
        std::cout << "# lambda=" << lambda << " energy_per_dimer=" << pcut::evaluate(result.energy_per_cell,lambda).real();
        if (particles) {
            const Eigen::SelfAdjointEigenSolver<pcut::Matrix> bands(result.bloch({momentum},lambda));
            if (bands.info() != Eigen::Success) throw std::runtime_error("Bloch eigensolver failed");
            std::cout << " k=" << momentum << " omega=" << bands.eigenvalues()[0];
        }
        std::cout << '\n';
        std::cerr << "Elapsed " << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() << " s\n";
    } catch (const std::exception& e) {
        std::cerr << "pcut_chain: " << e.what() << '\n';
        return 1;
    }
}
