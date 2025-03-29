#include <assert.h>

#include <numx/pde/fem.h>

int fem_pbc_slv(struct sim *sim, struct vec *q, struct fem_pbc_ops ops)
{
    (void)ops;

    assert(sim);
    assert(q);

    return 0;
}
