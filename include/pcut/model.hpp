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
    std::vector<int> charges; // basis zero is the unique local vacuum
    double vacuum_energy = 0.0;
    std::string name;
    void validate() const;
};
struct LocalTerm {
    std::vector<std::size_t> sites; // ordered tensor legs; first leg is least significant
    Matrix matrix;                // row = output, column = input
};

// H0 = sum(vacuum_energy) + gap * Q; V = sum(local terms), in physical units.
// Does not allocate the full tensor Hilbert space. Intermediate states are never
// truncated by their charge. User operators must include any statistics signs.
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
    [[nodiscard]] const std::vector<int>& changes() const noexcept { return changes_; }
    [[nodiscard]] SparseState apply(int change, const SparseState& state,
                                    std::size_t max_states = 1'000'000) const;
    [[nodiscard]] Matrix dense_hamiltonian(double lambda, State max_dimension = 4096) const;
private:
    struct Transition { State output; Complex value; };
    struct CompiledTerm {
        std::vector<std::size_t> sites;
        std::vector<State> local_stride;
        std::map<int, std::vector<std::vector<Transition>>> by_change;
    };
    std::vector<LocalSpace> spaces_;
    std::vector<State> strides_;
    std::vector<CompiledTerm> terms_;
    std::vector<int> changes_;
    State dimension_ = 1;
    double gap_;
};
[[nodiscard]] Complex evaluate(const Series& series, double lambda);
}
