#pragma once
#include <cstddef>
#include <string>
#include <vector>

namespace pcut { struct WhiteGraph; }

namespace pcut::detail {
// Simple undirected colored graph: physical vertices, edge occurrences,
// and ordered-leg ports have disjoint color namespaces. Channels color edges.
struct Incidence {
    std::vector<std::vector<int>> adjacent;
    std::vector<std::string> colors;
    int node(std::string color);
    void join(int a,int b);
};
// Inputs are validated by the caller; prepared signatures let construction
// reuse its exact interning pool. The convenience overload is used by tooling.
Incidence encode_incidence(const WhiteGraph& graph,
                           const std::vector<std::string>& labels,
                           const std::vector<std::string>& structures);
Incidence encode_incidence(const WhiteGraph& graph);
struct IncidenceLabel {
    std::vector<int> canonical_to_input;
    // Exact stabilizer indices, never the floating-point group-size statistics.
    std::vector<int> group_indices;
};
IncidenceLabel label_incidence(const Incidence& graph);
} // namespace pcut::detail
