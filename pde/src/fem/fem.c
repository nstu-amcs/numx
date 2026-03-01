#include <numx/pde/fem.h>
#include <numx/pde/sim.h>

#include "fem.h"

int fem_new(struct fem *fem)
{
    fem->slv.ops.tdd = TDD_I2S;
    fem->slv.ops.non.mod = NON_FPI;
    fem->slv.ops.non.ops.max = 200;
    fem->slv.ops.non.ops.err = 1e-10;
    fem->slv.ops.non.ops.rlx = false;
    fem->slv.ops.non.ops.itr.ctx = NULL;
    fem->slv.ops.non.ops.itr.run = NULL;
    fem->slv.ops.non.ops.dif = DIF_NUM;
    fem->slv.ops.non.ops.run.itr = 0;
    fem->slv.ops.non.ops.run.err = 0;
    fem->slv.ops.non.ops.run.rlx = 0;
    fem->slv.ops.iss.mod = ISS_BCG;
    fem->slv.ops.iss.ops.bcg.ops.max = 500;
    fem->slv.ops.iss.ops.bcg.ops.err = 1e-10;
    fem->slv.ops.iss.ops.bcg.ops.itr.ctx = NULL;
    fem->slv.ops.iss.ops.bcg.ops.itr.run = NULL;
    fem->slv.ops.iss.ops.bcg.ops.run.itr = 0;
    fem->slv.ops.iss.ops.bcg.ops.run.err = 0;
    fem->slv.ops.iss.ops.bcg.con.sm = NULL;
    fem->slv.exe = fem_exe;
    fem->slv.apx = NULL;
    fem->slv.itr.ctx = NULL;
    fem->slv.itr.run = NULL;
    fem->slv.run.wgt[0] = NULL;
    fem->slv.run.wgt[1] = NULL;
    fem->slv.run.wgt[2] = NULL;
    fem->slv.run.wgt[3] = NULL;
    fem->slv.run.bs = 1;
    fem->slv.run.ti = 0;
    fem->slv.run.tv = 0;
    fem->ops.bfs = FEM_BFS_LIN;

    return 0;
}

int fem_exe(struct sim *sim)
{
    switch (sim->mod) {
        case SIM_STD:
            return fem_std_exe(sim);
        case SIM_HMC:
            return fem_hmc_exe(sim);
    }

    return 0;
}
