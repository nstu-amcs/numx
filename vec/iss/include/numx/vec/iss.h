#ifndef NUMX_VEC_ISS_H
#define NUMX_VEC_ISS_H

#include <numx/vec/mtx.h>
#include <numx/vec/vec.h>
#include <stdx/cap.h>

enum iss_mod { ISS_JAC, ISS_RLX, ISS_BCG };

struct iss_itr {
  double err;
  int num;
};

struct iss_ops {
  /**
   * @brief Iteration callback closure.
   *
   * @param i current iteration
   * @param r current residual
   */
  struct vcap itr;

  /// @brief Target residual.
  double eps;

  /// @brief Maximum number of iterations.
  int max;
};

struct iss_jac_ops {
  struct iss_ops ops;

  /// @brief Relaxation parameter.
  double omg;
};

struct iss_rlx_ops {
  struct iss_ops ops;

  /// @brief Relaxation parameter.
  double omg;
};

struct iss_sor_ops {
  struct iss_ops ops;

  /// @brief Relaxation parameter.
  double omg;
};

struct iss_bcg_ops {
  struct iss_ops ops;

  /// @brief Preconditioner (maybe null).
  union con {
    struct smtx* sm;
  } con;
};

void iss_itr_cap(void* ctx, int n, ...);

int diss_jac_slv(struct dmtx* m, struct vec* x, struct vec* f, struct iss_jac_ops o);
int diss_rlx_slv(struct dmtx* m, struct vec* x, struct vec* f, struct iss_rlx_ops o);

int siss_bcg_slv(struct smtx* m, struct vec* x, struct vec* f, struct iss_bcg_ops o);

/// @brief Jacobi solver for systems of linear equations, Mx = f
#define iss_jac_slv(M, x, f, o) _Generic((M), \
  struct dmtx*: diss_jac_slv                  \
)(M, x, f, o)

/// @brief Relaxation (Gauss-Seidel) solver for systems of linear equations, Mx = f
#define iss_rlx_slv(M, x, f, o) _Generic((M), \
  struct dmtx*: diss_rlx_slv                  \
)(M, x, f, o)

/// @brief Successive Over-Relaxation solver for systems of linear equations, Mx = f
#define iss_sor_slv(M, x, f, o) _Generic((M), \
  struct dmtx*: diss_sor_slv                  \
)(M, x, f, o)

/// @brief Stabilized biconjugate gradient solver for systems of linear equations, Mx = f
#define iss_bcg_slv(M, x, f, o) _Generic((M), \
  struct smtx*: siss_bcg_slv                  \
)(M, x, f, o)

#endif  // NUMX_VEC_ISS_H
