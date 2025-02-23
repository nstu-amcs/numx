#include <errno.h>
#include <numx/pde/sse.h>
#include <numx/vec/iss.h>
#include <numx/vec/mtx.h>

int pde_sse_fdm_slv(struct obj* o, struct vec* x, struct fdm_sse_ops ops) {
  if (!o || !x) {
    errno = EINVAL;
    return -1;
  }

  struct dmtx m;
  struct vec f;

  if (mtx_new(&m, ((struct dmtx_pps){.n = o->v.len, .d = 5})))
    return -1;

  m.la[0] = 0;
  m.la[1] = 1;
  m.la[2] = 5;
  m.la[3] = -1;
  m.la[4] = -5;

  if (vec_new(&f, o->v.len))
    return -1;

  for (int i = 0; i < o->v.len; ++i) {
    struct vtx* v = o->v.dat[i];
    struct fdm_vtx_ctx* c = v->ctx;

    if (c->img) {
      m.ad[0][i] = 1;
      f.dat[i] = 1;

      continue;
    }

    switch (c->cnd.type) {
      case DIR:
        m.ad[0][i] = 1;
        f.dat[i] = c->cnd.pps.dir.tmp(v);

        continue;
      default:
        break;
    }

    double hl = v->x - ((struct vtx*)o->v.dat[i - 1])->x;
    double hr = ((struct vtx*)o->v.dat[i + 1])->x - v->x;
    double hu = ((struct vtx*)o->v.dat[i + ops.nx])->y - v->y;
    double hd = v->y - ((struct vtx*)o->v.dat[i - ops.nx])->y;

    m.ad[0][i] = 2 * ops.lam * (1 / (hl * hr) + 1 / (hd * hu)) + c->gam;

    m.ad[1][i] = -2 * ops.lam / (hr * (hr + hl));
    m.ad[2][i] = -2 * ops.lam / (hu * (hu + hd));
    m.ad[3][i] = -2 * ops.lam / (hl * (hr + hl));
    m.ad[4][i] = -2 * ops.lam / (hd * (hu + hd));

    f.dat[i] = c->ext;
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

  return 0;
}