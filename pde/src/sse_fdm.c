#include <errno.h>
#include <numx/pde/sse.h>
#include <numx/vec/iss.h>
#include <numx/vec/mtx.h>
#include <stdarg.h>
#include <stdio.h>

struct seg_ctx {
  struct cnd cnd;
};

struct qud_ctx {
  double lam;
  double gam;
};

int get_seg_ctx(void* ctx, int n, ...) {
  if (!ctx || n != 3) {
    errno = EINVAL;
    return -1;
  }

  struct pcut* dat = (struct pcut*)ctx;

  va_list arg;
  va_start(arg, n);

  const char* buf = va_arg(arg, const char*);
  void** sch = va_arg(arg, void**);

  va_end(arg);

  struct seg_ctx* sc = malloc(sizeof(struct seg_ctx));

  if (!sc)
    return -1;

  if (cnd_get(&sc->cnd, buf, dat)) {
    free(sc);
    return -1;
  }

  *sch = sc;

  return 0;
}

int get_qud_ctx(void* ctx, int n, ...) {
  (void)ctx;

  if (n != 3) {
    errno = EINVAL;
    return -1;
  }

  va_list arg;
  va_start(arg, n);

  const char* buf = va_arg(arg, const char*);
  void** qch = va_arg(arg, void**);

  va_end(arg);

  struct qud_ctx* qc = malloc(sizeof(struct qud_ctx));

  if (!qc)
    return -1;

  if (sscanf(buf, "%lf %lf", &qc->lam, &qc->gam) != 2) {
    free(qc);
    return -1;
  }

  *qch = qc;

  return 0;
}

enum seg_nrm { NRM_U, NRM_D, NRM_L, NRM_R };
enum vtx_sts { VTX_NON, VTX_DIR, VTX_NEU };

static enum seg_nrm seg_nrm(struct seg* seg, int nx) {
  int x1 = seg->vtx[0] % nx;
  int x2 = seg->vtx[1] % nx;
  int y1 = seg->vtx[0] / nx;
  int y2 = seg->vtx[1] / nx;

  if (y1 == y2)
    return x1 < x2 ? NRM_D : NRM_U;
  else
    return y1 < y2 ? NRM_R : NRM_L;
}

int pde_sse_fdm_slv(struct obj* obj, struct vec* x) {
  if (!obj || !x) {
    errno = EINVAL;
    return -1;
  }

  int r = 0;

  int nx = obj->ax.len;
  int ny = obj->ay.len;

  struct dmtx m;
  struct vec f;

  enum vtx_sts* sts = malloc(sizeof(enum vtx_sts) * nx * ny);

  if (!sts) {
    r = -1;
    goto end;
  }

  if ((r = mtx_new(&m, ((struct dmtx_pps){.n = nx * ny, .d = 5}))))
    goto end;

  if ((r = vec_new(&f, nx * ny)))
    goto end;

  for (int i = 0; i < nx * ny; ++i)
    sts[i] = VTX_NON;

  m.la[0] = 0;
  m.la[1] = 1;
  m.la[2] = nx;
  m.la[3] = -1;
  m.la[4] = -nx;

  for (int i = 0; i < nx * ny; ++i)
    m.ad[0][i] = 1;

  struct seg* seg = obj->seg.dat;

  for (int i = 0; i < obj->seg.len; ++i) {
    struct seg_ctx* ctx = seg[i].ctx;
    enum seg_nrm nrm = seg_nrm(&seg[i], nx);

    int step = 0;

    switch (nrm) {
      case NRM_L:
        step = nx;
        break;
      case NRM_R:
        step = -nx;
        break;
      case NRM_U:
        step = -1;
        break;
      case NRM_D:
        step = 1;
        break;
    }

    for (int i = seg[i].vtx[0]; i <= seg[i].vtx[1]; i += step) {
      struct vtx v = {.x = obj->ax.dat[i % nx], .y = obj->ay.dat[i / nx], .z = 0, .n = i, .ctx = 0};

      switch (sts[i]) {
        case VTX_NON:
          switch (ctx->cnd.type) {
            case DIR:
              m.ad[0][i] = 1;
              f.dat[i] = ctx->cnd.pps.dir.tmp(&v);
              sts[i] = VTX_DIR;

              break;
            case NEU:
              break;
            case ROB:
              break;
          }

          break;
        case VTX_DIR:
          break;
        case VTX_NEU:
          break;
      }
    }
  }

  // clang-format off
  if (iss_rlx_slv(&m, x, &f, ((struct iss_rlx_ops){
    .omg = 1,
    .ops = {
      .eps = 1e-10,
      .max = 1000,
    }
  }))) {
    mtx_cls(&m);
    vec_cls(&f);

    return -1;
  }
  // clang-format on

end:
  free(sts);

  mtx_cls(&m);
  vec_cls(&f);

  return r;
}