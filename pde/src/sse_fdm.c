#include <errno.h>
#include <numx/pde/sse.h>
#include <numx/vec/iss.h>
#include <numx/vec/mtx.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

typedef enum { NORM_U, NORM_D, NORM_L, NORM_R } norm;

typedef struct sctx {
  struct cnd cnd;
} sctx;

typedef struct qctx {
  double lam;
  double gam;

  double (*f)(struct vtx*);
} qctx;

int sctx_get(void* ctx, int n, ...) {
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

  sctx* sc = malloc(sizeof(sctx));

  if (!sc)
    return -1;

  if (cnd_get(&sc->cnd, buf, dat)) {
    free(sc);
    return -1;
  }

  *sch = sc;

  return 0;
}

int qctx_get(void* ctx, int n, ...) {
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

  qctx* qc = malloc(sizeof(qctx));

  if (!qc)
    return -1;

  if (sscanf(buf, "%lf %lf", &qc->lam, &qc->gam) != 2) {
    free(qc);
    return -1;
  }

  *qch = qc;

  return 0;
}

static norm seg_norm(struct seg* s, int nx) {
  int x1 = s->vtx[0] % nx;
  int x2 = s->vtx[1] % nx;
  int y1 = s->vtx[0] / nx;
  int y2 = s->vtx[1] / nx;

  if (y1 == y2)
    return x1 < x2 ? NORM_D : NORM_U;
  else
    return y1 < y2 ? NORM_R : NORM_L;
}

static int seg_evo(struct obj* obj, struct dmtx* m, struct vec* f, bool* status) {
  int nx = obj->ax.len;
  int ns = obj->seg.len;

  struct seg* ss = (struct seg*)obj->seg.dat;

  for (int i = 0; i < ns; ++i) {
    sctx* ctx = ss[i].ctx;
    norm norm = seg_norm(&ss[i], nx);

    int i0 = 0;
    int i1 = 0;
    int st = 0;

    switch (norm) {
      case NORM_L:
        i0 = ss[i].vtx[1];
        i1 = ss[i].vtx[0];
        st = nx;
        break;
      case NORM_R:
        i0 = ss[i].vtx[0];
        i1 = ss[i].vtx[1];
        st = nx;
        break;
      case NORM_U:
        i0 = ss[i].vtx[1];
        i1 = ss[i].vtx[0];
        st = 1;
        break;
      case NORM_D:
        i0 = ss[i].vtx[0];
        i1 = ss[i].vtx[1];
        st = 1;
        break;
    }

    for (int i = i0; i <= i1; i += st) {
      if (status[i])
        continue;

      struct vtx v = {.x = obj->ax.dat[i % nx], .y = obj->ay.dat[i / nx], .z = 0, .n = i, .ctx = 0};

      switch (ctx->cnd.type) {
        case DIR:
          m->ad[0][i] = 1;
          f->dat[i] = ctx->cnd.pps.dir.tmp(&v);

          break;
        case NEU:
          f->dat[i] = ctx->cnd.pps.neu.tta(&v);

          double hi = 0;

          switch (norm) {
            case NORM_L:
              hi = obj->ax.dat[i + 1] - obj->ax.dat[i];

              m->ad[0][i] = 1 / hi;
              m->ad[1][i] = -1 / hi;

              break;
            case NORM_R:
              hi = obj->ax.dat[i] - obj->ax.dat[i - 1];

              m->ad[0][i] = 1 / hi;
              m->ad[3][i] = -1 / hi;

              break;
            case NORM_U:
              hi = obj->ay.dat[i] - obj->ay.dat[i - nx];

              m->ad[0][i] = 1 / hi;
              m->ad[4][i] = -1 / hi;

              break;
            case NORM_D:
              hi = obj->ay.dat[i + nx] - obj->ay.dat[i];

              m->ad[0][i] = 1 / hi;
              m->ad[2][i] = -1 / hi;

              break;
          }

          break;
        case ROB:
          break;
      }

      status[i] = true;
    }
  }

  return 0;
}

static int qud_evo(struct obj* obj, struct dmtx* m, struct vec* f, bool* status) {
  int nx = obj->ax.len;
  int nq = obj->qud.len;

  double* ax = obj->ax.dat;
  double* ay = obj->ay.dat;

  for (int i = 0; i < nq; ++i) {
    int a = obj->qud.dat[i].vtx[0];
    int b = obj->qud.dat[i].vtx[1];
    int c = obj->qud.dat[i].vtx[2];

    int dx = nx - (b - a) - 1;

    int x0 = a % nx;
    int x1 = b % nx;
    int y0 = a / nx;
    int y1 = c / nx;

    qctx* ctx = obj->qud.dat[i].ctx;

    for (int i = y0, n = a; i <= y1; ++i, n += dx) {
      for (int j = x0; j <= x1; ++j, ++n) {
        if (status[n])
          continue;

        struct vtx v = {.x = ax[j], .y = ay[i], .z = 0, .n = n, .ctx = 0};

        f->dat[n] = ctx->f(&v);

        double hu = ay[i + 1] - ay[i];
        double hr = ax[j + 1] - ax[j];
        double hd = ay[i] - ay[i - 1];
        double hl = ax[j] - ax[j - 1];

        m->ad[0][n] = 2 * ctx->lam * (1 / (hl * hr) + 1 / (hd * hu)) + ctx->gam;
        m->ad[1][n] = -2 * ctx->lam / (hr * (hr + hl));
        m->ad[2][n] = -2 * ctx->lam / (hu * (hu + hd));
        m->ad[3][n] = -2 * ctx->lam / (hl * (hr + hl));
        m->ad[4][n] = -2 * ctx->lam / (hd * (hu + hd));
      }
    }
  }

  return 0;
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

  bool* status = malloc(sizeof(bool) * nx * ny);

  if (!status) {
    r = -1;
    goto end;
  }

  if ((r = mtx_new(&m, ((struct dmtx_pps){.n = nx * ny, .d = 5}))))
    goto end;

  if ((r = vec_new(&f, nx * ny)))
    goto end;

  memset(status, 0, sizeof(bool) * nx * ny);

  m.la[0] = 0;
  m.la[1] = 1;
  m.la[2] = nx;
  m.la[3] = -1;
  m.la[4] = -nx;

  for (int i = 0; i < nx * ny; ++i)
    m.ad[0][i] = 1;

  if ((r = seg_evo(obj, &m, &f, status)))
    goto end;

  if ((r = qud_evo(obj, &m, &f, status)))
    goto end;

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
  free(status);

  mtx_cls(&m);
  vec_cls(&f);

  return r;
}