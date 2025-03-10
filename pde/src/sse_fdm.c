#include <errno.h>
#include <numx/pde/sse.h>
#include <numx/vec/iss.h>
#include <numx/vec/mtx.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

enum stat { S_NON, S_INN, S_DIR, S_NEU, S_ROB };

static enum norm seg_norm(struct seg* s, int nx) {
  int x1 = s->vtx[0] % nx;
  int x2 = s->vtx[1] % nx;
  int y1 = s->vtx[0] / nx;
  int y2 = s->vtx[1] / nx;

  if (y1 == y2)
    return x1 < x2 ? NORM_D : NORM_U;
  else
    return y1 < y2 ? NORM_R : NORM_L;
}

static int seg_evo(struct obj* obj, struct dmtx* m, struct vec* f, enum stat* status) {
  int nx = obj->ax.len;
  int ns = obj->seg.len;

  struct seg* ss = (struct seg*)obj->seg.dat;

  for (int i = 0; i < ns; ++i) {
    struct fctx* ctx = ss[i].ctx;
    enum norm norm = seg_norm(&ss[i], nx);

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
      int xi = i % nx;
      int yi = i / nx;

      struct vtx v = {.x = obj->ax.dat[xi], .y = obj->ay.dat[yi], .z = 0, .n = i, .ctx = 0};

      double hi = 0;
      int bi = 0;

      switch (ctx->cnd.type) {
        case DIR:
          m->ad[i][0] = 1;
          f->dat[i] = ctx->cnd.pps.dir.tmp(&v);

          status[i] = S_DIR;

          break;
        case NEU:
          if (status[i] == S_DIR || status[i] == S_ROB)
            continue;

          f->dat[i] += ctx->cnd.pps.neu.tta(&v);

          hi = 0;
          bi = 0;

          switch (norm) {
            case NORM_L:
              hi = obj->ax.dat[xi + 1] - obj->ax.dat[xi];
              bi = 1;

              break;
            case NORM_R:
              hi = obj->ax.dat[xi] - obj->ax.dat[xi - 1];
              bi = 3;

              break;
            case NORM_U:
              hi = obj->ay.dat[yi] - obj->ay.dat[yi - 1];
              bi = 4;

              break;
            case NORM_D:
              hi = obj->ay.dat[yi + 1] - obj->ay.dat[yi];
              bi = 2;

              break;
          }

          m->ad[i][0] += 1 / hi;
          m->ad[i][bi] += -1 / hi;

          if (status[i] == S_NEU) {
            f->dat[i] /= 2;
            m->ad[i][0] /= 2;
            m->ad[i][bi] /= 2;
          }

          status[i] = S_NEU;

          break;
        case ROB:
          if (status[i] == S_DIR || status[i] == S_NEU)
            continue;

          f->dat[i] += ctx->cnd.pps.rob.bet * ctx->cnd.pps.rob.tmp(&v);

          hi = 0;
          bi = 0;

          switch (norm) {
            case NORM_L:
              hi = obj->ax.dat[xi + 1] - obj->ax.dat[xi];
              bi = 1;

              break;
            case NORM_R:
              hi = obj->ax.dat[xi] - obj->ax.dat[xi - 1];
              bi = 3;

              break;
            case NORM_U:
              hi = obj->ay.dat[yi] - obj->ay.dat[yi - 1];
              bi = 4;

              break;
            case NORM_D:
              hi = obj->ay.dat[yi + 1] - obj->ay.dat[yi];
              bi = 2;

              break;
          }

          m->ad[i][0] += 1 / hi + ctx->cnd.pps.rob.bet;
          m->ad[i][bi] += -1 / hi;

          if (status[i] == S_ROB) {
            f->dat[i] /= 2;
            m->ad[i][0] /= 2;
            m->ad[i][bi] /= 2;
          }

          status[i] = S_ROB;

          break;
      }
    }
  }

  return 0;
}

static int qud_evo(struct obj* obj, struct dmtx* m, struct vec* f, enum stat* status) {
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

    struct ectx* ctx = obj->qud.dat[i].ctx;

    for (int i = y0, n = a; i <= y1; ++i, n += dx) {
      if (status[n] && status[n + 1]) {
        n += b - a + 1;
        continue;
      }

      for (int j = x0; j <= x1; ++j, ++n) {
        if (status[n])
          continue;

        struct vtx v = {.x = ax[j], .y = ay[i], .z = 0, .n = n, .ctx = 0};

        f->dat[n] = ctx->ext(&v);

        double hu = ay[i + 1] - ay[i];
        double hr = ax[j + 1] - ax[j];
        double hd = ay[i] - ay[i - 1];
        double hl = ax[j] - ax[j - 1];

        m->ad[n][0] = 2 * ctx->lam * (1 / (hl * hr) + 1 / (hd * hu)) + ctx->gam;
        m->ad[n][1] = -2 * ctx->lam / (hr * (hr + hl));
        m->ad[n][2] = -2 * ctx->lam / (hu * (hu + hd));
        m->ad[n][3] = -2 * ctx->lam / (hl * (hr + hl));
        m->ad[n][4] = -2 * ctx->lam / (hd * (hu + hd));

        status[n] = S_INN;
      }
    }
  }

  return 0;
}

struct itr {
  int k;
  double r;
};

void cback(void* ctx, int n, ...) {
  va_list arg;
  va_start(arg, n);

  ((struct itr*)ctx)->k = va_arg(arg, int);
  ((struct itr*)ctx)->r = va_arg(arg, double);

  va_end(arg);
}

int pde_sse_fdm_slv(struct obj* obj, struct vec* x, struct sse_fdm_ops ops) {
  if (!obj || !x) {
    errno = EINVAL;
    return -1;
  }

  int r = 0;

  int nx = obj->ax.len;
  int ny = obj->ay.len;

  struct dmtx m;
  struct vec f;

  enum stat* status = malloc(sizeof(enum stat) * nx * ny);

  if (!status) {
    r = -1;
    goto end;
  }

  memset(status, 0, sizeof(bool) * nx * ny);

  if ((r = mtx_new(&m, ((struct dmtx_pps){.n = nx * ny, .d = 5}))))
    goto end;

  if ((r = vec_new(&f, nx * ny)))
    goto end;

  if ((r = vec_new(x, nx * ny)))
    goto end;

  m.la[0] = 0;
  m.la[1] = 1;
  m.la[2] = nx;
  m.la[3] = -1;
  m.la[4] = -nx;

  if ((r = seg_evo(obj, &m, &f, status)))
    goto end;

  if ((r = qud_evo(obj, &m, &f, status)))
    goto end;

  for (int i = 0; i < nx * ny; ++i)
    if (!status[i])
      m.ad[i][0] = 1;

  switch (ops.ops.iss.mod) {
    case ISS_JAC:
      if ((r = iss_jac_slv(&m, x, &f, ops.ops.iss.ops.jac)))
        goto end;

      break;
    case ISS_RLX:
      if ((r = iss_rlx_slv(&m, x, &f, ops.ops.iss.ops.rlx)))
        goto end;

      break;
    default:
      r = -1;
      errno = ENOTSUP;
      goto end;
  }

end:
  free(status);

  mtx_cls(&m);
  vec_cls(&f);

  return r;
}