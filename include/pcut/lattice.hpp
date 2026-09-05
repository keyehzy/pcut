#pragma once
#include <pcut/model.hpp>
#include <compare>
#include <cstddef>
#include <vector>

namespace pcut {
using Coordinate = std::vector<int>;
struct Site {
    Coordinate cell;
    std::size_t basis = 0;
    auto operator<=>(const Site&) const = default;
};
struct OperatorChannel {
    Matrix matrix; // fixed Hermitian operator, in physical energy units
    bool fermionic = false;
};
struct CoupledChannel {
    OperatorChannel op;
    double coupling = 1; // dimensionless ratio; all channels share one lambda
};
struct Interaction {
    std::vector<Site> legs; // ordered offsets and unit-cell basis labels
    std::vector<CoupledChannel> channels;
};
struct PeriodicLattice {
    unsigned dimension = 1;
    std::vector<LocalSpace> cell;
    std::vector<Interaction> interactions; // distinct physical hyperedge templates
    double gap = 1.0;
    void validate() const;
    [[nodiscard]] std::vector<std::vector<double>> couplings() const;
};
struct Edge {
    std::size_t type;
    Coordinate origin;
    auto operator<=>(const Edge&) const = default;
};
using Cluster = std::vector<Edge>;
struct NormalizedCluster {
    Cluster cluster;
    Coordinate shift; // original = normalized + shift
};
[[nodiscard]] NormalizedCluster normalize(Cluster cluster);
[[nodiscard]] std::vector<Site> vertices(const PeriodicLattice& lattice, const Cluster& cluster);
[[nodiscard]] bool connected(const PeriodicLattice& lattice, const Cluster& cluster);
[[nodiscard]] ClusterModel cluster_model(const PeriodicLattice& lattice, const Cluster& cluster);
}
