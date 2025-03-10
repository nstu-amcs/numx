#ifndef NUMX_PDE_SSE_H
#define NUMX_PDE_SSE_H

#include <numx/pde/cnd.h>
#include <numx/pde/geo.h>
#include <numx/vec/iss.h>
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

struct sse_ops {
  struct {
    enum iss_mod mod;
    union {
      struct iss_jac_ops jac;
      struct iss_rlx_ops rlx;
      struct iss_bcg_ops bcg;
    } ops;
  } iss;
};

struct sse_fdm_ops {
  struct sse_ops ops;
};

struct sse_fem_ops {
  enum { FEM_LIN } mod;

  struct sse_ops ops;
};

/// @brief Finite Difference Method boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fdm_slv(struct obj* obj, struct vec* x, struct sse_fdm_ops ops);

/// @brief Finite Element Method boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fem_slv(struct obj* obj, struct vec* x, struct sse_fem_ops ops);

#endif  // NUMX_PDE_SSE_H