#include "fem.h"

int fem_ell_lin_non_slv(struct sim *sim, struct vec *wgt)
{
    int r = 0;

    int    max = sim->fem->ops.non.max;
    double err = sim->fem->ops.non.err;
    double cur = 1;

    struct vec t;

    double tn = 0;
    double bn = 0;

    if ((r = vec_new(&t, sim->fem->vec.n)))
        goto end;

    for (int i = 0; i < max && cur > err; ++i) {
        if ((r = fem_ell_lin_std_slv(sim, wgt)))
            goto end;

        mtx_vmlt(&sim->fem->mtx, wgt, &t);
        vec_cmb(&t, &sim->fem->vec, &t, -1);
        vec_nrm(&t, &tn);
        vec_nrm(&sim->fem->vec, &bn);

        cur = tn / bn;
    }

end:
    vec_cls(&t);

    return r;
}
