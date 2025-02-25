#include <errno.h>
#include <numx/pde/sse.h>
#include <numx/vec/iss.h>
#include <numx/vec/mtx.h>

int pde_sse_fdm_slv(struct obj* o, struct vec* x) {
  if (!o || !x) {
    errno = EINVAL;
    return -1;
  }

  struct fdm_obj_ctx* oc = o->ctx;

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
    struct fdm_vtx_ctx* vc = v->ctx;

    if (vc->img) {
      m.ad[0][i] = 1;
      f.dat[i] = 1;

      continue;
    }

    switch (vc->cnd.type) {
      case DIR:
        m.ad[0][i] = 1;
        f.dat[i] = vc->cnd.pps.dir.tmp(v);

        continue;
      default:
        break;
    }

    double hr = ((struct vtx*)o->v.dat[i + 1])->x - v->x;
    double hu = ((struct vtx*)o->v.dat[i + o->pps.nx])->y - v->y;
    double hl = v->x - ((struct vtx*)o->v.dat[i - 1])->x;
    double hd = v->y - ((struct vtx*)o->v.dat[i - o->pps.nx])->y;

    m.ad[0][i] = 2 * oc->lam * (1 / (hl * hr) + 1 / (hd * hu)) + vc->gam;

    m.ad[1][i] = -2 * oc->lam / (hr * (hr + hl));
    m.ad[2][i] = -2 * oc->lam / (hu * (hu + hd));
    m.ad[3][i] = -2 * oc->lam / (hl * (hr + hl));
    m.ad[4][i] = -2 * oc->lam / (hd * (hu + hd));

    f.dat[i] = vc->ext;
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