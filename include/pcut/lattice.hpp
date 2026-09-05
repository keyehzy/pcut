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
struct Interaction {
    std::vector<Site> legs; // ordered offsets and unit-cell basis labels
    Matrix matrix;
    bool fermionic = false; // even operator in ordered-leg Fock basis
};
struct PeriodicLattice {
    unsigned dimension = 1;
    std::vector<LocalSpace> cell;
    std::vector<Interaction> interactions; // distinct colored hyperedge templates
    double gap = 1.0;
    void validate() const;
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
struct Subcluster {
    std::size_t index;
    std::vector<std::size_t> vertex_map; // normalized child's vertices -> parent's vertices
};
struct ClusterEntry {
    Cluster edges;
    std::vector<Site> sites;
    std::vector<Subcluster> subclusters; // ALL proper connected embedded subsets, including multiplicity
};
struct EnumerationOptions {
    std::size_t max_clusters = 100'000;
    std::size_t max_subclusters = 5'000'000;
};
// Translation classes of actual colored embedded bond animals (no white graphs).
// Each entry has one embedding per unit cell; rotations/reflections stay distinct.
class ClusterCatalog {
public:
    ClusterCatalog(PeriodicLattice lattice, unsigned max_edges, EnumerationOptions options = {});
    [[nodiscard]] const PeriodicLattice& lattice() const noexcept { return lattice_; }
    [[nodiscard]] const std::vector<ClusterEntry>& entries() const noexcept { return entries_; }
    [[nodiscard]] unsigned max_edges() const noexcept { return max_edges_; }
private:
    PeriodicLattice lattice_;
    unsigned max_edges_;
    std::vector<ClusterEntry> entries_;
};
}
