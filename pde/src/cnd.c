#include <assert.h>
#include <errno.h>
#include <stdio.h>

#include <numx/pde/cnd.h>

cut_gen(fun_cut, fun, STDX_PUB);

static const char *skip(const char *buf)
{
    while (*buf == ' ' || *buf == '\t' || *buf == '|')
        buf++;

    return buf;
}

int cnd_get(struct cnd *cnd, const char *buf, struct fun_cut *dat)
{
    assert(cnd);
    assert(buf);
    assert(dat);

    buf = skip(buf);

    char t = *buf;
    int  f = 0;
    int  n = 0;

    buf = skip(buf + 1);

    if (sscanf(buf, "%d%n", &f, &n) != 1)
        return -1;

    if (f < 0 || f >= dat->len)
        return -1;

    switch (t) {
        case 'D':
            cnd->mod = CND_DIR;
            cnd->pps.dir.tmp = dat->dat[f];

            break;
        case 'N':
            cnd->mod = CND_NEU;
            cnd->pps.neu.tta = dat->dat[f];

            break;
        case 'R':
            cnd->mod = CND_ROB;
            cnd->pps.rob.tmp = dat->dat[f];

            buf = skip(buf + n);

            if (sscanf(buf, "%d", &f) != 1)
                return -1;

            if (f < 0 || f >= dat->len)
                return -1;

            cnd->pps.rob.bet = dat->dat[f];

            break;
        default:
            errno = EINVAL;
            return -1;
    }

    return 0;
}
