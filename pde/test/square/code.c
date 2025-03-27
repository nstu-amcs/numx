#include <numx/pde/sim.h>

double diffusion(struct sim *sim, int vtx, double tgt)
{
    (void)sim;
    (void)vtx;
    (void)tgt;

    return 5;
}

double source(struct sim *sim, int vtx, double tgt)
{
    (void)sim;
    (void)vtx;
    (void)tgt;

    return 0;
}
