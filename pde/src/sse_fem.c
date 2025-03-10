#include <assert.h>
#include <numx/pde/sse.h>
#include <numx/vec/iss.h>
#include <numx/vec/mtx.h>

static int qud_evo();
static int hxd_evo();

int pde_sse_fem_slv(struct obj* obj, struct vec* x, struct sse_fem_ops ops) {
  assert(obj);
  assert(x);

  int r = 0;
  int n = obj->ax.len * obj->ay.len * obj->az.len;

  struct smtx m;
  struct vec f;

  if (mtx_new(&m, (struct smtx_pps){.n = n}))
    goto end;

  if (vec_new(&f, n))
    goto end;

  if (hxd_evo())
    goto end;

  if (qud_evo())
    goto end;

  switch (ops.ops.iss.mod) {
    case ISS_BCG:
      if ((r = iss_bcg_slv(&m, x, &f, ops.ops.iss.ops.bcg)))
        goto end;

      break;

    default:
      r = -1;
      errno = ENOTSUP;
      goto end;
  }

end:
  mtx_cls(&m);
  vec_cls(&f);

  return 0;
}

static int qud_evo() { return 0; }

static int hxd_evo() { return 0; }