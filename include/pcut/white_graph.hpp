#pragma once
#include <pcut/lattice.hpp>
#include <pcut/effective.hpp>
#include <functional>
#include <mutex>
#include <optional>

namespace pcut {
struct GraphEdge {
    std::vector<std::size_t> legs;
    std::vector<OperatorChannel> channels;
};
struct WhiteGraph {
    std::vector<LocalSpace> spaces; // name identifies species; basis order is significant
    std::vector<GraphEdge> edges; // parallel edges are distinct occurrences
    [[nodiscard]] ClusterModel model(double gap) const;
    [[nodiscard]] std::size_t variables() const;
};
struct GraphMap {
    std::vector<std::size_t> vertices; // source -> target
    std::vector<std::size_t> edges;
    std::vector<std::size_t> channels; // flattened edge/channel source -> target
};
struct CanonicalGraph {
    WhiteGraph graph;
    std::string key; // exact structural serialization, never a hash-only identity
    GraphMap map; // canonical -> supplied graph
    std::size_t vertex_automorphisms = 0;
};
[[nodiscard]] CanonicalGraph canonicalize(const WhiteGraph& graph);
struct GraphSubcluster {
    std::size_t graph;
    GraphMap map; // child's canonical graph -> parent canonical graph
};
struct GraphEntry {
    CanonicalGraph canonical;
    std::vector<GraphSubcluster> subclusters; // all connected edge subsets, including multiplicity
};
struct EmbeddedSubcluster {
    std::size_t index;
    GraphMap map; // physical child -> physical parent
};
struct GraphEmbedding {
    std::size_t graph;
    Cluster edges; // normalized physical edge set; one occurrence per translation class
    std::vector<Site> sites; // sorted physical sites
    GraphMap map; // canonical -> physical sorted vertices, edges, flattened channels
    std::vector<EmbeddedSubcluster> subclusters;
};
// Own this object across sweeps. Callbacks must be pure, isomorphism-covariant,
// cluster additive and return homogeneous corrections with degree n at order n.
// Captured evaluation context must remain immutable. Copies share identity.
class ScalarEvaluator {
public:
    using Function=std::function<SymbolicSeries(const WhiteGraph&, double, unsigned)>;
    explicit ScalarEvaluator(Function function);
    [[nodiscard]] const std::shared_ptr<const Function>& function() const { return function_; }
private:
    std::shared_ptr<const Function> function_;
};
class GraphCache {
public:
    [[nodiscard]] std::shared_ptr<const SymbolicBlock> block(const CanonicalGraph& graph,
        double gap, const EffectiveOperator& effective, const std::vector<State>& basis, bool linked = false);
    [[nodiscard]] std::shared_ptr<const SymbolicSeries> scalar(const CanonicalGraph& graph,
        double gap, unsigned order, const ScalarEvaluator& evaluator);
    [[nodiscard]] std::size_t evaluations() const;
    [[nodiscard]] std::size_t hits() const;
private:
    friend class WhiteGraphExpansion;
    void record_hit();
    mutable std::mutex mutex_;
    std::map<std::string,std::shared_ptr<const SymbolicBlock>> blocks_;
    std::map<std::string,std::shared_ptr<const SymbolicSeries>> scalars_;
    // Retain callback ownership, so identity addresses cannot be recycled.
    std::vector<std::shared_ptr<const ScalarEvaluator::Function>> evaluators_;
    std::size_t evaluations_=0, hits_=0;
};
// Canonical abstract graphs collected from exhaustive infinite embedding witnesses.
// Geometry and numerical couplings are absent from canonical graph identity.
class WhiteGraphExpansion {
public:
    WhiteGraphExpansion(PeriodicLattice lattice, unsigned max_edges,
                        std::shared_ptr<GraphCache> cache = std::make_shared<GraphCache>());
    [[nodiscard]] const PeriodicLattice& lattice() const noexcept { return lattice_; }
    [[nodiscard]] unsigned max_edges() const noexcept { return max_edges_; }
    [[nodiscard]] const std::vector<GraphEntry>& graphs() const noexcept;
    [[nodiscard]] const std::vector<GraphEmbedding>& embeddings() const noexcept;
    // Ratios in physical interaction/channel order. Shares the immutable plan.
    [[nodiscard]] WhiteGraphExpansion bind(const std::vector<std::vector<double>>& couplings) const;
    [[nodiscard]] const std::vector<double>& couplings(const GraphEmbedding& embedding) const;
    // Compiled unbound operators in physical site order; ratios are not substituted.
    [[nodiscard]] const ClusterModel& structural_model(const GraphEmbedding& embedding) const;
    [[nodiscard]] const std::shared_ptr<GraphCache>& cache() const noexcept { return cache_; }
    // Pass an embedding reference from this expansion. External basis uses sorted physical sites. Both state and formal variable
    // relabelings are applied, including the graded Fock permutation.
    // linked=true projects onto full-edge monomial support; vacuum subtraction
    // is still required for a one-particle kernel.
    [[nodiscard]] std::vector<Matrix> block(const GraphEmbedding& embedding,
        const EffectiveOperator& effective, const std::vector<State>& basis, bool linked = false) const;
private:
    PeriodicLattice lattice_;
    unsigned max_edges_;
    std::shared_ptr<GraphCache> cache_;
    struct Plan;
    std::shared_ptr<const Plan> plan_;
    std::vector<std::vector<double>> couplings_;
};
}
