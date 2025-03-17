#include <assert.h>
#include <errno.h>
#include <numx/pde/sse.h>
#include <stdarg.h>

static const char *skip(const char *buf)
{
    while (*buf == ' ' || *buf == '\t' || *buf == '|')
        buf++;

    return buf;
}

int fctx_get(void *ctx, int n, ...)
{
    assert(ctx);
    assert(n == 2);

    struct fun_cut *dat = (struct fun_cut *)ctx;

    va_list arg;
    va_start(arg, n);

    const char *buf = va_arg(arg, const char *);
    void      **fch = va_arg(arg, void **);

    va_end(arg);

    struct fctx *fc = malloc(sizeof(struct fctx));

    if (!fc)
        return -1;

    if (cnd_get(&fc->cnd, buf, dat)) {
        free(fc);
        return -1;
    }

    *fch = fc;

    return 0;
}

int ectx_get(void *ctx, int n, ...)
{
    assert(ctx);
    assert(n == 2);

    struct fun_cut *dat = (struct fun_cut *)ctx;

    va_list arg;
    va_start(arg, n);

    const char *buf = va_arg(arg, const char *);
    void      **ech = va_arg(arg, void **);

    va_end(arg);

    struct ectx *ec = malloc(sizeof(struct ectx));

    if (!ec)
        return -1;

    int f = 0;
    int p = 0;

    buf = skip(buf);

    if (sscanf(buf, "%d%n", &f, &p) != 1) {
        free(ec);
        return -1;
    }

    if (f < 0 || f >= dat->len) {
        free(ec);
        errno = EINVAL;
        return -1;
    }

    ec->lam = dat->dat[f];
    buf = skip(buf + p);

    if (sscanf(buf, "%d%n", &f, &p) != 1) {
        free(ec);
        return -1;
    }

    if (f < 0 || f >= dat->len) {
        free(ec);
        errno = EINVAL;
        return -1;
    }

    ec->gam = dat->dat[f];
    buf = skip(buf + p);

    if (sscanf(buf, "%d", &f) != 1) {
        free(ec);
        return -1;
    }

    if (f < 0 || f >= dat->len) {
        free(ec);
        errno = EINVAL;
        return -1;
    }

    ec->ext = dat->dat[f];
    *ech = ec;

    return 0;
}
