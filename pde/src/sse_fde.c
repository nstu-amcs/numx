#include <errno.h>
#include <numx/pde/sse.h>

int pde_sse_fde_slv(struct obj* o, struct vec* x) {
  if (!o || !x) {
    errno = EINVAL;
    return -1;
  }

  return 0;
}