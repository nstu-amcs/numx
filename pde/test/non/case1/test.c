#include <math.h>
#include <numx/pde/sim.h>
#include <stdio.h>

#define N 3

double target(struct vtx *v) { return v->x + v->y + v->z; }

void on_slv(void *ctx, struct sim *sim) {
  struct vtx *vtx = sim->msh->vtx.dat;

  double err = 0;

  for (int h = 0; h < sim->msh->hxd.len; ++h) {
    struct hxd *hxd = &sim->msh->hxd.dat[h];

    int v0 = hxd->vtx[0];
    int v7 = hxd->vtx[7];

    double x0 = vtx[v0].x;
    double y0 = vtx[v0].y;
    double z0 = vtx[v0].z;

    double x1 = vtx[v7].x;
    double y1 = vtx[v7].y;
    double z1 = vtx[v7].z;

    double xs = (x1 - x0) / N;
    double ys = (y1 - y0) / N;
    double zs = (z1 - z0) / N;

    struct vtx v = {x0, y0, z0};

    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        for (int k = 0; k < N; ++k) {
          double apx =
              sim->slv->apx.run(sim->slv->apx.ctx, sim, sim->rt.wgt, h, &v);
          double tgt = target(&v);

          err += fabs(tgt - apx);
          v.z += zs;
        }

        v.y += ys;
      }

      v.x += xs;
    }
  }

  *((double *)ctx) = err;
}

int main(int argc, char **argv) {
  (void)argc;
  (void)argv;

  struct sim sim;

  sim_new(&sim);
  sim_imp_elm(&sim, "elm/case.sif");

  struct fem *fem = (struct fem *)sim.slv;
  double err = 0;

  fem->slv.itr.ctx = &err;
  fem->slv.itr.run = on_slv;

  sim_run(&sim);

  printf("Simulation:\n\tError: %.7e\n\tNonlinear:\n\t\tIterations: "
         "%d\n\t\tError: %.7e",
         err, fem->ops.non.rt.itr, fem->ops.non.rt.err);

  sim_cls(&sim);

  return 0;
}
