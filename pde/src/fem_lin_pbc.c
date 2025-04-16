#include <assert.h>

#include <numx/pde/sim.h>
#include <numx/vec/vec.h>

#include "fem_lin.h"

static int fem_lin_pbc_i1s(struct sim *sim, struct vec *wgt);

int fem_lin_pbc_slv(struct sim *sim) {
  assert(sim);

  int r = 0;
  int n = sim->msh->vtx.len;

  struct fem *fem = (struct fem *)sim->slv;
  struct vec w1;
  struct vec w2;

  if ((r = vec_new(&w1, n)))
    goto end;

  if ((r = vec_new(&w2, n)))
    goto end;

  if ((r = fem_lin_mtx_lam_asm(sim, &fem->prv.ell)))
    goto end;

  if ((r = fem_lin_mtx_gam_asm(sim, &fem->prv.ell)))
    goto end;

  if ((r = fem_lin_mtx_rob_asm(sim, &fem->prv.ell)))
    goto end;

  if ((r = fem_lin_mtx_sig_asm(sim, &fem->prv.pbc)))
    goto end;

  int num = sim->ops.tdd.num;
  int hop = sim->ops.tdd.hop;

  for (int i = 0; i < n; ++i)
    fem->prv.ell.dr[i] += fem->prv.pbc.dr[i] / hop;

  for (int i = 0; i < fem->prv.pbc.pps.z; ++i) {
    fem->prv.ell.lr[i] += fem->prv.pbc.lr[i] / hop;
    fem->prv.ell.ur[i] += fem->prv.pbc.ur[i] / hop;
  }

  if (fem_lin_mtx_dir_asm(sim, &fem->prv.ell))
    goto end;

  if (fem_lin_pbc_i1s(sim, &w1))
    goto end;

  sim->rt.ti = 0;
  sim->rt.tv = 0;
  sim->rt.wgt = &w1;

  if (fem->slv.itr.run)
    fem->slv.itr.run(fem->slv.itr.ctx, sim);

  if (sim->ops.exp.put.run)
    sim->ops.exp.put.run(sim->ops.exp.put.ctx, sim);

  for (int i = 1, t = hop; i <= num; i += 1, t += hop) {
    sim->rt.ti = i;
    sim->rt.tv = t;
    sim->rt.wgt = &w2;

    if ((r = mtx_vmlt(&fem->prv.pbc, &w1, &fem->prv.vec)))
      goto end;

    for (int j = 0; j < n; ++j)
      fem->prv.vec.dat[j] /= hop;

    if ((r = fem_lin_vec_src_asm(sim, &fem->prv.vec)))
      goto end;

    if ((r = fem_lin_vec_neu_asm(sim, &fem->prv.vec)))
      goto end;

    if ((r = fem_lin_vec_rob_asm(sim, &fem->prv.vec)))
      goto end;

    if ((r = fem_lin_vec_dir_asm(sim, &fem->prv.vec)))
      goto end;

    switch (fem->ops.iss.mod) {
    case ISS_BCG:
      if (iss_bcg_slv(&fem->prv.ell, &w2, &fem->prv.vec, &fem->ops.iss.ops.bcg))
        return -1;

      break;
    default:
      errno = ENOTSUP;
      return -1;
    }

    if (fem->slv.itr.run)
      fem->slv.itr.run(fem->slv.itr.ctx, sim);

    if (sim->ops.exp.put.run)
      sim->ops.exp.put.run(sim->ops.exp.put.ctx, sim);

    w1.dat = w2.dat;
  }

end:
  vec_cls(&w1);
  vec_cls(&w2);

  return r;
}

static int fem_lin_pbc_i1s(struct sim *sim, struct vec *wgt) {
  assert(sim);
  assert(wgt);

  for (int i = 0; i < sim->msh->hxd.len; ++i) {
    struct hxd *hxd = &sim->msh->hxd.dat[i];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct cnd_ini *ini = &sim->cnd_ini.dat[obj->ini];

    if (ini->tgt.type == VAL_FUN)
      for (int k = 0; k < 8; ++k)
        wgt->dat[hxd->vtx[k]] = ini->tgt.as.fun(sim, hxd->vtx[k]);
    else
      for (int k = 0; k < 8; ++k)
        wgt->dat[hxd->vtx[k]] = ini->tgt.as.num;
  }

  return 0;
}
