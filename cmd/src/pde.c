#include <stdio.h>

#include <numx/cmd.h>
#include <numx/pde/fem.h>
#include <numx/pde/sim.h>

static void man()
{
    printf("usage: numx pde {filename.sif}\n");
}

int pde(int argc, char **argv)
{
    if (argc < 3) {
        man();
        return -1;
    }

    struct sim sim;

    if (sim_new(&sim, argv[2])) {
        perror("simulation");
        return -1;
    }

    struct vec res;

    if (vec_new(&res, sim.msh->vtx.len)) {
        perror("memory");
        return -1;
    }

    if (fem_slv(&sim, &res)) {
        perror("solver");
        return -1;
    }

    return 0;
}
