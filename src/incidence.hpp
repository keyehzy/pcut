#pragma once
#include <cstddef>
#include <string>
#include <vector>

namespace pcut::detail {
// Simple undirected colored graph: physical vertices, edge occurrences,
// ordered-leg ports, and channel occurrences have disjoint color namespaces.
struct Incidence {
    std::vector<std::vector<int>> adjacent;
    std::vector<std::string> colors;
    int node(std::string color);
    void join(int a,int b);
};
struct IncidenceLabel {
    std::vector<int> canonical_to_input;
    // Exact stabilizer indices, never the floating-point group-size statistics.
    std::vector<int> group_indices;
};
IncidenceLabel label_incidence(const Incidence& graph);
} // namespace pcut::detail
