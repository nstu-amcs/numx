#ifndef NUMX_VEC_ISS_H
#define NUMX_VEC_ISS_H

#include <numx/vec/mtx.h>
#include <numx/vec/vec.h>
#include <stdx/cap.h>

struct iss_opt {
  /**
   * @brief Iteration callback closure
   *
   * @param i current iteration
   * @param r current residual
   */
  struct vcap itr;

  /// @brief Target residual
  double eps;

  /// @brief Maximum number of iterations
  int max;
};

struct iss_jac_opt {
  struct iss_opt opt;

  /// @brief Relaxation parameter
  double omg;
};

struct iss_rlx_opt {
  struct iss_opt opt;

  /// @brief Relaxation parameter
  double omg;
};

struct iss_bcg_opt {
  struct iss_opt opt;

  /// @brief Preconditioner (maybe null)
  union con {
    struct smtx* sm;
  } con;
};

int diss_jac_slv(struct dmtx* m, struct vec* x, struct vec* f, struct iss_jac_opt o);
int diss_rlx_slv(struct dmtx* m, struct vec* x, struct vec* f, struct iss_rlx_opt o);

int siss_bcg_slv(struct smtx* m, struct vec* x, struct vec* f, struct iss_bcg_opt o);

/// @brief Jacobi solver for systems of linear equations, Mx = f
#define iss_jac_slv(M, x, f, o) _Generic((M), \
  struct dmtx*: diss_jac_slv                  \
)(M, x, f, o)

/// @brief Relaxation (Gauss-Seidel) solver for systems of linear equations, Mx = f
#define iss_rlx_slv(M, x, f, o) _Generic((M), \
  struct dmtx*: diss_rlx_slv                  \
)(M, x, f, o)

/// @brief Stabilized biconjugate gradient solver for systems of linear equations, Mx = f
#define iss_bcg_slv(M, x, f, o) _Generic((M), \
  struct smtx*: siss_bcg_slv                  \
)(M, x, f, o)

#endif  // NUMX_VEC_ISS_H
