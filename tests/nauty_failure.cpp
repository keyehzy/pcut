#include <pcut/pcut.hpp>
#include <cstddef>
#include <new>
#include <iostream>
#include <stdexcept>
extern "C" void pcut_nauty_fail_after(std::ptrdiff_t);
extern "C" bool pcut_nauty_allocations_empty();
int main() {
    pcut::WhiteGraph graph;
    graph.spaces.resize(9,pcut::models::ising_chain().cell[0]);
    for (std::size_t v=1;v<9;++v)
        graph.edges.push_back({{0,v},{{pcut::Matrix::Identity(4,4),false}}});
    const auto expected=pcut::canonicalize(graph).key;
    bool completed=false;
    std::size_t failures=0;
    for (std::ptrdiff_t position=0;!completed;++position) {
        pcut_nauty_fail_after(position);
        try {
            if (pcut::canonicalize(graph).key!=expected)
                throw std::runtime_error("incorrect canonical key during fault injection");
            completed=true;
        }
        catch (const std::bad_alloc&) { ++failures; }
        pcut_nauty_fail_after(-1);
        if (!pcut_nauty_allocations_empty()) throw std::runtime_error("leaked upstream allocation");
        if (pcut::canonicalize(graph).key!=expected) throw std::runtime_error("failed recovery");
    }
    if (!failures) throw std::runtime_error("no upstream allocation failures injected");
    std::cout << "Injected " << failures << " upstream allocation failures\n";
}
