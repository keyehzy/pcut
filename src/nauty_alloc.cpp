#include <cstdlib>
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
namespace {
// Track upstream-owned allocations without allocating tracking storage. This
// also recovers a tcnode if its second allocation fails before nauty links it.
struct alignas(std::max_align_t) Allocation {
    Allocation *previous,*next;
    std::size_t bytes;
};
thread_local Allocation* allocations=nullptr;
#ifdef PCUT_NAUTY_ALLOC_TESTING
thread_local std::ptrdiff_t failure_after=-1;
#endif
}
#ifdef PCUT_NAUTY_ALLOC_TESTING
extern "C" void pcut_nauty_fail_after(std::ptrdiff_t count) { failure_after=count; }
extern "C" bool pcut_nauty_allocations_empty() { return allocations==nullptr; }
#endif
extern "C" void* pcut_nauty_allocate(std::size_t count,std::size_t size) {
#ifdef PCUT_NAUTY_ALLOC_TESTING
    if (failure_after==0) throw std::bad_alloc();
    if (failure_after>0) --failure_after;
#endif
    if (size && count>(std::numeric_limits<std::size_t>::max()-sizeof(Allocation))/size)
        throw std::length_error("nauty allocation size overflow");
    auto* block=static_cast<Allocation*>(std::malloc(sizeof(Allocation)+count*size));
    if (!block) throw std::bad_alloc();
    block->previous=nullptr; block->next=allocations; block->bytes=count*size;
    if (allocations) allocations->previous=block;
    allocations=block; return block+1;
}
extern "C" void pcut_nauty_free(void* p) {
    if (!p) return;
    auto* block=static_cast<Allocation*>(p)-1;
    if (block->previous) block->previous->next=block->next;
    else allocations=block->next;
    if (block->next) block->next->previous=block->previous;
    std::free(block);
}
extern "C" void* pcut_nauty_reallocate(void* old,std::size_t size) {
    auto* replacement=pcut_nauty_allocate(1,size);
    if (old) std::memcpy(replacement,old,std::min(size,(static_cast<Allocation*>(old)-1)->bytes));
    pcut_nauty_free(old); return replacement;
}
extern "C" void pcut_nauty_release_unowned() {
    while (allocations) pcut_nauty_free(allocations+1);
}
