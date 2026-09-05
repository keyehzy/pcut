#include <pcut/lattice.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace pcut {
namespace {
Coordinate shift(const Coordinate& a, const Coordinate& b, int sign = 1) {
    if (a.size() != b.size()) throw std::invalid_argument("coordinate dimension mismatch");
    Coordinate c(a.size());
    for (std::size_t d = 0; d < a.size(); ++d) {
        const auto v = static_cast<long long>(a[d]) + static_cast<long long>(sign)*b[d];
        if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max())
            throw std::overflow_error("lattice coordinate overflow");
        c[d] = static_cast<int>(v);
    }
    return c;
}
std::vector<Site> edge_vertices(const PeriodicLattice& l, const Edge& e) {
    if (e.type >= l.interactions.size() || e.origin.size() != l.dimension)
        throw std::invalid_argument("invalid embedded edge");
    auto sites = l.interactions[e.type].legs;
    for (auto& s : sites) s.cell = shift(s.cell, e.origin);
    return sites;
}
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
        // Reuse all operator and Hilbert-space validation.
        const ClusterModel checked(std::move(spaces), {{std::move(legs), interaction.matrix}}, gap);
        (void)checked;
    }
}
NormalizedCluster normalize(Cluster cluster) {
    if (cluster.empty()) throw std::invalid_argument("cannot normalize an empty edge cluster");
    Coordinate origin = cluster.front().origin;
    for (const auto& e : cluster) origin = std::min(origin, e.origin);
    for (auto& e : cluster) e.origin = shift(e.origin, origin, -1);
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
        t.matrix = lattice.interactions.at(e.type).matrix;
        terms.push_back(std::move(t));
    }
    return {std::move(spaces), std::move(terms), lattice.gap};
}
ClusterCatalog::ClusterCatalog(PeriodicLattice lattice, unsigned max_edges, EnumerationOptions options)
    : lattice_(std::move(lattice)), max_edges_(max_edges) {
    lattice_.validate();
    if (max_edges > 64) throw std::invalid_argument("cluster order exceeds 64");
    std::set<Cluster> level;
    if (max_edges > 0) for (std::size_t t = 0; t < lattice_.interactions.size(); ++t)
        level.insert({{t, Coordinate(lattice_.dimension,0)}});
    std::map<Cluster, std::size_t> indices;
    std::size_t subcluster_count = 0;
    for (unsigned size = 1; size <= max_edges && !level.empty(); ++size) {
        if (level.size() > options.max_clusters - entries_.size()) throw std::length_error("cluster budget exceeded");
        for (const auto& c : level) {
            ClusterEntry entry{c, vertices(lattice_,c), {}};
            // Recursive deletion generates all connected subsets without a 2^n bit-mask limit.
            std::set<Cluster> visited;
            std::vector<Cluster> pending{c};
            for (std::size_t p = 0; p < pending.size(); ++p) {
                const Cluster parent = pending[p];
                if (parent.size() == 1) continue;
                for (std::size_t erase = 0; erase < parent.size(); ++erase) {
                    Cluster child = parent;
                    child.erase(child.begin()+static_cast<std::ptrdiff_t>(erase));
                    if (!connected(lattice_,child) || !visited.insert(child).second) continue;
                    if (++subcluster_count > options.max_subclusters) throw std::length_error("subcluster budget exceeded");
                    pending.push_back(child);
                    const auto normalized = normalize(child);
                    const auto index = indices.at(normalized.cluster);
                    Subcluster sub{index,{}};
                    for (auto s : entries_[index].sites) {
                        s.cell = shift(s.cell, normalized.shift);
                        const auto it = std::lower_bound(entry.sites.begin(), entry.sites.end(), s);
                        if (it == entry.sites.end() || *it != s) throw std::logic_error("subcluster embedding mismatch");
                        sub.vertex_map.push_back(static_cast<std::size_t>(it-entry.sites.begin()));
                    }
                    entry.subclusters.push_back(std::move(sub));
                }
            }
            indices.emplace(c, entries_.size());
            entries_.push_back(std::move(entry));
        }
        if (size == max_edges) break;
        std::set<Cluster> next;
        for (const auto& c : level) for (const auto& site : vertices(lattice_,c)) {
            for (std::size_t t = 0; t < lattice_.interactions.size(); ++t)
                for (const auto& leg : lattice_.interactions[t].legs) if (leg.basis == site.basis) {
                    const Edge e{t, shift(site.cell, leg.cell, -1)};
                    if (std::binary_search(c.begin(),c.end(),e)) continue;
                    auto enlarged = c;
                    enlarged.push_back(e);
                    next.insert(normalize(std::move(enlarged)).cluster);
                    if (next.size() > options.max_clusters - entries_.size()) throw std::length_error("cluster budget exceeded");
                }
        }
        level = std::move(next);
    }
}
}
