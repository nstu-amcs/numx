#include <string.h>

#include <numx/pde/sim.h>

int main() {
  struct sim sim;

  sim_new(&sim);

  strcpy(sim.ops.exp.pfx, "u1");
  strcpy(sim.ops.usr.pfx, "code");

  sim_imp_elm(&sim, "pde/test/pbc/case4/elm/case.sif");

  sim_run(&sim);
  sim_cls(&sim);

  return 0;
}
