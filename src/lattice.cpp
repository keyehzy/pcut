#include <pcut/lattice.hpp>
#include "detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace pcut {
namespace {
std::vector<Site> edge_vertices(const PeriodicLattice& l, const Edge& e) {
    if (e.type >= l.interactions.size() || e.origin.size() != l.dimension)
        throw std::invalid_argument("invalid embedded edge");
    auto sites = l.interactions[e.type].legs;
    for (auto& s : sites) s.cell = detail::translate(s.cell, e.origin);
    return sites;
}
}
std::vector<std::vector<double>> PeriodicLattice::couplings() const {
    std::vector<std::vector<double>> result;
    for (const auto& interaction : interactions) {
        result.emplace_back();
        for (const auto& channel : interaction.channels) result.back().push_back(channel.coupling);
    }
    return result;
}
void PeriodicLattice::validate() const {
    if (dimension == 0 || dimension > 16 || cell.empty() || !std::isfinite(gap) || gap <= 0)
        throw std::invalid_argument("invalid lattice dimension, unit cell, or gap");
    for (const auto& s : cell) s.validate();
    for (const auto& interaction : interactions) {
        std::vector<LocalSpace> spaces;
        std::vector<std::size_t> legs;
        std::set<Site> unique;
        for (const auto& site : interaction.legs) {
            if (site.basis >= cell.size() || site.cell.size() != dimension || !unique.insert(site).second)
                throw std::invalid_argument("invalid or repeated interaction site");
            legs.push_back(spaces.size());
            spaces.push_back(cell[site.basis]);
        }
        if (interaction.channels.empty()) throw std::invalid_argument("interaction needs operator channels");
        for (const auto& channel : interaction.channels) {
            if (!std::isfinite(channel.coupling)) throw std::invalid_argument("nonfinite coupling");
            const ClusterModel checked(spaces, {{legs, channel.op.matrix, channel.op.fermionic}}, gap);
            (void)checked;
        }
    }
}
NormalizedCluster normalize(Cluster cluster) {
    if (cluster.empty()) throw std::invalid_argument("cannot normalize an empty edge cluster");
    Coordinate origin = cluster.front().origin;
    for (const auto& e : cluster) origin = std::min(origin, e.origin);
    for (auto& e : cluster) e.origin = detail::translate(e.origin, origin, -1);
    std::sort(cluster.begin(), cluster.end());
    if (std::adjacent_find(cluster.begin(), cluster.end()) != cluster.end())
        throw std::invalid_argument("duplicate embedded interaction");
    return {std::move(cluster), std::move(origin)};
}
std::vector<Site> vertices(const PeriodicLattice& lattice, const Cluster& cluster) {
    std::set<Site> sites;
    for (const auto& e : cluster) for (auto s : edge_vertices(lattice,e)) sites.insert(std::move(s));
    return {sites.begin(), sites.end()};
}
bool connected(const PeriodicLattice& lattice, const Cluster& cluster) {
    if (cluster.empty()) return false;
    std::set<Site> visited;
    std::vector<bool> seen(cluster.size(), false);
    seen[0] = true;
    for (auto s : edge_vertices(lattice, cluster[0])) visited.insert(std::move(s));
    bool changed = true;
    while (changed) {
        changed = false;
        for (std::size_t i = 0; i < cluster.size(); ++i) if (!seen[i]) {
            const auto sites = edge_vertices(lattice, cluster[i]);
            if (std::any_of(sites.begin(), sites.end(), [&](const Site& s) { return visited.contains(s); })) {
                seen[i] = true; changed = true;
                visited.insert(sites.begin(), sites.end());
            }
        }
    }
    return std::all_of(seen.begin(), seen.end(), [](bool x) { return x; });
}
ClusterModel cluster_model(const PeriodicLattice& lattice, const Cluster& cluster) {
    const auto sites = vertices(lattice, cluster);
    std::vector<LocalSpace> spaces;
    for (const auto& s : sites) spaces.push_back(lattice.cell.at(s.basis));
    std::vector<LocalTerm> terms;
    for (const auto& e : cluster) {
        LocalTerm t;
        for (const auto& s : edge_vertices(lattice,e))
            t.sites.push_back(static_cast<std::size_t>(std::lower_bound(sites.begin(), sites.end(),s)-sites.begin()));
        for (const auto& channel : lattice.interactions.at(e.type).channels) {
            t.matrix=channel.coupling*channel.op.matrix;
            t.fermionic=channel.op.fermionic;
            terms.push_back(t);
        }
    }
    return {std::move(spaces), std::move(terms), lattice.gap};
}
}
