#include <errno.h>
#include <numx/pde/sse.h>
#include <stdarg.h>

int fctx_get(void* ctx, int n, ...) {
  if (!ctx || n != 2) {
    errno = EINVAL;
    return -1;
  }

  struct ocut* dat = (struct ocut*)ctx;

  va_list arg;
  va_start(arg, n);

  const char* buf = va_arg(arg, const char*);
  void** fch = va_arg(arg, void**);

  va_end(arg);

  struct fctx* fc = malloc(sizeof(struct fctx));

  if (!fc)
    return -1;

  if (cnd_get(&fc->cnd, buf, dat)) {
    free(fc);
    return -1;
  }

  *fch = fc;

  return 0;
}

int ectx_get(void* ctx, int n, ...) {
  if (!ctx || n != 2) {
    errno = EINVAL;
    return -1;
  }

  struct ocut* dat = (struct ocut*)ctx;

  va_list arg;
  va_start(arg, n);

  const char* buf = va_arg(arg, const char*);
  void** ech = va_arg(arg, void**);

  va_end(arg);

  struct ectx* ec = malloc(sizeof(struct ectx));

  if (!ec)
    return -1;

  int f = 0;

  if (sscanf(buf, "%lf %lf %d", &ec->lam, &ec->gam, &f) != 3) {
    free(ec);
    return -1;
  }

  if (f >= dat->len) {
    free(ec);

    errno = EINVAL;
    return -1;
  }

  ec->ext = dat->dat[f];
  *ech = ec;

  return 0;
}