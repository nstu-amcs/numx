#include <numx/pde/sim.h>

double diffusion(struct sim *sim, int vtx, double tgt)
{
    (void)sim;
    (void)vtx;

    return 2 * tgt;
}

double reaction(struct sim *sim, int vtx, double tgt)
{
    (void)sim;
    (void)vtx;

    return 2 * tgt;
}
