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

    if (sim_new(&sim)) {
        perror("initialization");
        return -1;
    }

    if (sim_imp_elm(&sim, argv[2])) {
        perror("import");
        return -1;
    }

    if (sim_run(&sim)) {
        perror("execution");
        return -1;
    }

    sim_cls(&sim);

    return 0;
}
