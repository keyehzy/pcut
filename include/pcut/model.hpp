#pragma once
#include <Eigen/Core>
#include <complex>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace pcut {
using Complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using State = std::uint64_t;
using SparseState = std::unordered_map<State, Complex>;
using Series = std::vector<Complex>; // coefficient of lambda^n at index n

struct LocalSpace {
    std::vector<int> charges; // nonnegative integer Q; basis zero has Q=0
    double vacuum_energy = 0.0;
    std::string name;
    std::vector<int> particles; // optional conserved electron number per local state
    std::vector<unsigned> parity; // optional Z2 grading, 0 even / 1 odd
    void validate() const;
    void require_product_vacuum() const;
};
struct LocalTerm {
    std::vector<std::size_t> sites; // ordered tensor legs; first leg is least significant
    Matrix matrix;                // row = output, column = input
    bool fermionic = false;       // even operator in ordered-leg Fock basis
};

// H0 = sum(vacuum_energy) + gap * Q; V = sum(local terms), in physical units.
// Does not allocate the full tensor Hilbert space. Intermediate states are never
// truncated by their charge. Fermionic terms receive graded embedding signs.
class ClusterModel {
public:
    ClusterModel(std::vector<LocalSpace> spaces, std::vector<LocalTerm> terms,
                 double gap = 1.0, double zero_tolerance = 0.0);
    [[nodiscard]] const std::vector<LocalSpace>& spaces() const noexcept { return spaces_; }
    [[nodiscard]] std::size_t sites() const noexcept { return spaces_.size(); }
    [[nodiscard]] State dimension() const noexcept { return dimension_; }
    [[nodiscard]] double gap() const noexcept { return gap_; }
    [[nodiscard]] double vacuum_energy() const noexcept;
    [[nodiscard]] State encode(const std::vector<unsigned>& local) const;
    [[nodiscard]] std::vector<unsigned> decode(State state) const;
    [[nodiscard]] int charge(State state) const;
    [[nodiscard]] int particle_number(State state) const;
    void require_product_vacuum() const;
    // Reject ambiguous tensor/graded mixtures before operator-valued linking.
    void require_operator_grading() const;
    [[nodiscard]] bool conserves_particles() const noexcept { return conserves_particles_; }
    [[nodiscard]] bool fermionic() const noexcept { return fermionic_; }
    [[nodiscard]] const std::vector<int>& changes() const noexcept { return changes_; }
    [[nodiscard]] SparseState apply(int change, const SparseState& state,
                                    std::size_t max_states = 1'000'000) const;
    [[nodiscard]] Matrix dense_hamiltonian(double lambda, State max_dimension = 4096) const;
private:
    struct Transition { State output; Complex value; };
    struct CompiledTerm {
        std::vector<std::size_t> sites;
        std::vector<State> local_stride;
        bool fermionic = false;
        std::vector<std::pair<std::size_t, std::size_t>> inversions;
        std::map<int, std::vector<std::vector<Transition>>> by_change;
    };
    std::vector<LocalSpace> spaces_;
    std::vector<State> strides_;
    std::vector<CompiledTerm> terms_;
    std::vector<int> changes_;
    State dimension_ = 1;
    double gap_;
    bool conserves_particles_ = true;
    bool fermionic_ = false;
    bool operator_grading_valid_ = true;
};
// Sign of gathering ordered support sites to the front, with all other sites
// following in their original order. Used on BOTH input and output states.
[[nodiscard]] int graded_gather_sign(const std::vector<LocalSpace>& spaces,
                                     const std::vector<unsigned>& local,
                                     const std::vector<std::size_t>& support);
[[nodiscard]] Complex evaluate(const Series& series, double lambda);
}
