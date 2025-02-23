#ifndef NUMX_PDE_SSE_H
#define NUMX_PDE_SSE_H

#include <numx/pde/cnd.h>
#include <numx/pde/geo.h>
#include <numx/vec/vec.h>

/// @brief Vertex context for Finite Difference Method solver.
struct fdm_vtx_ctx {
  /// @brief Wether vertex imaginary or not.
  bool img;

  /// @brief External influence.
  double ext;

  /// @brief Gamma coefficient.
  double gam;

  /// @brief Boundary condition.
  struct cnd cnd;
};

/// @brief Options for solid-state equations.
struct sse_ops {
  /// @brief Diffusion coefficient.
  double lam;
};

/// @brief Finite Difference Method boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fdm_slv(struct obj* o, struct vec* x, struct sse_ops ops);

/// @brief Finite Element Method boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fem_slv(struct obj* o, struct vec* x);

#endif  // NUMX_PDE_SSE_H