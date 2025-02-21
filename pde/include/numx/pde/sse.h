#ifndef NUMX_PDE_SSE_H
#define NUMX_PDE_SSE_H

#include <numx/pde/cnd.h>
#include <numx/pde/geo.h>
#include <numx/vec/vec.h>

/// @brief Vertex context for Finite Difference Method solver.
struct vtx_fde_ctx {
  /// @brief Wether vertex imaginary or not.
  bool img;

  /// @brief Boundary condition.
  struct cnd cnd;
};

/// @brief Finite Difference Method boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fde_slv(struct obj* o, struct vec* x);

/// @brief Finite Element Method boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fem_slv(struct obj* o, struct vec* x);

#endif  // NUMX_PDE_SSE_H