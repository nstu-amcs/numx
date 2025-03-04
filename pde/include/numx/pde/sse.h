#ifndef NUMX_PDE_SSE_H
#define NUMX_PDE_SSE_H

#include <numx/pde/cnd.h>
#include <numx/pde/geo.h>
#include <numx/vec/vec.h>
#include <stdx/cut.h>

/// @brief Face's context.
struct fctx {
  struct cnd cnd;
};

int fctx_get(void* ctx, int n, ...);

/// @brief Element's context.
struct ectx {
  double lam;
  double gam;

  double (*ext)(struct vtx* v);
};

int ectx_get(void* ctx, int n, ...);

/// @brief Finite Difference Method boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fdm_slv(struct obj* obj, struct vec* x);

/// @brief Finite Element Method boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fem_slv(struct obj* obj, struct vec* x);

#endif  // NUMX_PDE_SSE_H