#pragma once
#include "incidence.hpp"
#include <nausparse.h>

namespace pcut::detail {
// Own all sparse input/output arrays in C++; views never own their pointers.
struct SparseStorage {
    std::vector<std::size_t> starts;
    std::vector<int> degrees, neighbors;
    sparsegraph view();
};
struct PreparedIncidence {
    SparseStorage graph;
    std::vector<int> lab, partition;
    explicit PreparedIncidence(const Incidence& input);
};
// Also runs during exception unwinding. Traces-specific scratch is tooling-only.
struct NautyScratch {
    NautyScratch()=default;
    NautyScratch(const NautyScratch&)=delete;
    NautyScratch& operator=(const NautyScratch&)=delete;
    ~NautyScratch();
};
} // namespace pcut::detail
