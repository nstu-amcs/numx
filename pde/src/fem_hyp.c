#include <assert.h>
#include <numx/pde/fem.h>

int fem_hyp_slv(struct sim *sim, struct vec *q, struct fem_hyp_ops ops)
{
    (void)ops;

    assert(sim);
    assert(q);

    return 0;
}
