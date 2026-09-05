// Force-included only into upstream C translation units. Allocate before freeing
// old scratch, so a C++ allocation exception leaves freedyn() safe to call.
#include <nauty.h>
#include <stdlib.h>
extern void *pcut_nauty_allocate(size_t, size_t);
extern void *pcut_nauty_reallocate(void *, size_t);
extern void pcut_nauty_free(void *);
#undef ALLOCS
#undef REALLOCS
#undef FREES
#define ALLOCS(x,y) pcut_nauty_allocate((size_t)(x),(size_t)(y))
#define REALLOCS(p,x) pcut_nauty_reallocate(p,(size_t)(x))
#define FREES(p) pcut_nauty_free(p)
#undef DYNALLOC1
#undef DYNALLOC2
#undef DYNREALLOC
#define DYNALLOC1(type,name,name_sz,sz,msg) \
 if ((size_t)(sz) > (name_sz)) { \
   type *pcut_new = (type*)pcut_nauty_allocate((size_t)(sz),sizeof(type)); \
   FREES(name); name=pcut_new; name_sz=(size_t)(sz); }
#define DYNALLOC2(type,name,name_sz,sz1,sz2,msg) \
 DYNALLOC1(type,name,name_sz,(size_t)(sz1)*(size_t)(sz2),msg)
#define DYNREALLOC(type,name,name_sz,sz,msg) \
 if ((size_t)(sz) > (name_sz)) { \
   name=(type*)pcut_nauty_reallocate(name,(size_t)(sz)*sizeof(type)); \
   name_sz=(size_t)(sz); }
