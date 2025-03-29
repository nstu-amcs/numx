#include <assert.h>
#include <cgnslib.h>

#include <numx/pde/msh.h>

int msh_imp_gns(struct msh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);

    errno = ENOTSUP;

    return -1;
}

int msh_exp_gns(struct msh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);

    return 0;
}
