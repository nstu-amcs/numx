#include <assert.h>
#include <stdarg.h>

#include <numx/com/cmp.h>

int iasc(void *ctx, int n, ...)
{
    (void)ctx;

    assert(n == 2);

    va_list arg;
    va_start(arg, n);

    int a = va_arg(arg, int);
    int b = va_arg(arg, int);

    va_end(arg);

    if (a < b)
        return 1;

    if (b < a)
        return -1;

    return 0;
}

int idsc(void *ctx, int n, ...)
{
    (void)ctx;

    assert(n == 2);

    va_list arg;
    va_start(arg, n);

    int a = va_arg(arg, int);
    int b = va_arg(arg, int);

    va_end(arg);

    if (b < a)
        return 1;

    if (a < b)
        return -1;

    return 0;
}
