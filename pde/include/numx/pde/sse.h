#ifndef NUMX_PDE_SSE_H
#define NUMX_PDE_SSE_H

#include <numx/pde/cnd.h>
#include <numx/pde/geo.h>
#include <numx/vec/vec.h>

/// @brief Vertex context for Finite Difference Method solver.
struct fdm_vtx_ctx {
  /// @brief Wether vertex imaginary or not.
  bool img;

  /// @brief Boundary condition.
  struct cnd cnd;

  double ext;
  double gam;
};

struct fdm_sse_ops {
  double lam;

  int nx;
  int ny;
};

struct fem_sse_ops {};

/// @brief Finite Difference Method boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fdm_slv(struct obj* o, struct vec* x, struct fdm_sse_ops ops);

/// @brief Finite Element Method boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fem_slv(struct obj* o, struct vec* x, struct fem_sse_ops ops);

#endif  // NUMX_PDE_SSE_H