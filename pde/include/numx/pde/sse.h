#ifndef NUMX_PDE_SSE_H
#define NUMX_PDE_SSE_H

#include <numx/pde/cnd.h>
#include <numx/pde/geo.h>
#include <numx/vec/vec.h>

int pde_sse_fem_slv(struct obj* o, struct vec* x);
int pde_sse_fde_slv(struct obj* o, struct vec* x);

#endif  // NUMX_PDE_SSE_H