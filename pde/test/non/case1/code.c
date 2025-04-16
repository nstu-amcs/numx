#include <numx/pde/sim.h>

double source(struct sim *sim, int vtx) {
  (void)sim;
  (void)vtx;

  struct vtx *v = &sim->msh->vtx.dat[vtx];

  double x = v->x;
  double y = v->y;
  double z = v->z;

  return -30 + 0.4 * (x + y + z);
}

double lam(struct sim *sim, int vtx) {
  (void)sim;
  (void)vtx;

  struct vtx *v = &sim->msh->vtx.dat[vtx];

  double x = v->x;
  double y = v->y;
  double z = v->z;

  return 10 * (x + y + z) + 1;
}

double neu1(struct sim *sim, int vtx) {
  (void)sim;
  (void)vtx;

  struct vtx *v = &sim->msh->vtx.dat[vtx];

  double x = v->x;
  double y = v->y;
  double z = v->z;

  return -10 * (x + y + z) - 1;
}

double neu2(struct sim *sim, int vtx) {
  (void)sim;
  (void)vtx;

  struct vtx *v = &sim->msh->vtx.dat[vtx];

  double x = v->x;
  double y = v->y;
  double z = v->z;

  return 10 * (x + y + z) + 1;
}

double neu3(struct sim *sim, int vtx) {
  (void)sim;
  (void)vtx;

  struct vtx *v = &sim->msh->vtx.dat[vtx];

  double x = v->x;
  double y = v->y;
  double z = v->z;

  return 10 * (x + y + z) + 1;
}

double rob1(struct sim *sim, int vtx) {
  (void)sim;
  (void)vtx;

  struct vtx *v = &sim->msh->vtx.dat[vtx];

  double x = v->x;
  double y = v->y;
  double z = v->z;

  return 20.1 + x + 2 * y + 2 * z;
}

double rob2(struct sim *sim, int vtx) {
  (void)sim;
  (void)vtx;

  struct vtx *v = &sim->msh->vtx.dat[vtx];

  double y = v->y;

  return -y - 0.1;
}

double rob3(struct sim *sim, int vtx) {
  (void)sim;
  (void)vtx;

  struct vtx *v = &sim->msh->vtx.dat[vtx];

  double z = v->z;

  return -z - 0.1;
}
