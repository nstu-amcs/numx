#include <assert.h>
#include <dlfcn.h>
#include <stdio.h>

#include <numx/pde/sim.h>
#include <numx/pde/val.h>

static int get(const char *buf, char *key, char *val);

static int hdr_get(FILE *f, struct sim *sim);
static int sim_get(FILE *f, struct sim *sim);
static int obj_get(FILE *f, struct sim *sim);
static int mat_get(FILE *f, struct sim *sim);
static int ext_get(FILE *f, struct sim *sim);
static int ini_get(FILE *f, struct sim *sim);
static int bnd_get(FILE *f, struct sim *sim);

int sim_new(struct sim *sim, const char *sif)
{
    assert(sim);
    assert(sif);

    FILE *f = fopen(sif, "r");

    if (!f)
        return -1;

    char buf[128];

    while (fgets(buf, sizeof(buf), f)) {
        switch (*buf) {
            case 'H':
                if (hdr_get(f, sim))
                    return -1;

                break;
            case 'S':
                if (sim_get(f, sim))
                    return -1;

                break;
            case 'B':
                if (!strncmp(buf, "Body Force", 10)) {
                    if (ext_get(f, sim))
                        return -1;

                    continue;
                }

                if (!strncmp(buf, "Boundary Condition", 18)) {
                    if (bnd_get(f, sim))
                        return -1;

                    continue;
                }

                if (!strncmp(buf, "Body", 4)) {
                    if (obj_get(f, sim))
                        return -1;

                    continue;
                }

                break;
            case 'M':
                if (mat_get(f, sim))
                    return -1;

                break;
            case 'I':
                if (ini_get(f, sim))
                    return -1;

                break;
        }
    }

    fclose(f);

    return 0;
}

static const char *path_get(const char *src, char *dst);

static int hdr_get(FILE *f, struct sim *sim)
{
    char buf[128];
    char path[64];

    while (fgets(buf, sizeof(buf), f)) {
        const char *cur = buf;

        while (*cur == ' ')
            ++cur;

        switch (*cur) {
            case 'E':
                return 0;
            case 'M':
                cur = path_get(cur, sim->pps.msh.dir);
                cur = path_get(cur, sim->pps.msh.pfx);

                break;
            case 'I':
                cur = path_get(cur, path);
                sim->pps.usr = dlopen(path, 1);

                break;
            case 'R':
                cur = path_get(cur, sim->pps.exp.dir);

                break;
        }
    }

    return 0;
}

static const char *path_get(const char *src, char *dst)
{
    while (*src != '"')
        ++src;

    src += 1;

    for (int i = 0; *src != '"'; ++i, ++src)
        dst[i] = *src;

    return src + 1;
}

static int get(const char *buf, char *key, char *val)
{
    while (*buf == ' ')
        ++buf;

    const char *beg = buf;
    const char *act = buf;

    while (*buf != '=') {
        if (*buf != ' ')
            act = buf;

        ++buf;
    }

    strncpy(key, beg, act - beg);

    while (*buf == ' ')
        ++buf;

    beg = buf;
    act = buf;

    while (*buf != '\n' && *buf != '\0') {
        if (*buf != ' ')
            act = buf;

        ++buf;
    }

    strncpy(val, beg, act - beg);

    return 0;
}

static int sim_get(FILE *f, struct sim *sim)
{
    return 0;
}

static int obj_get(FILE *f, struct sim *sim);
static int mat_get(FILE *f, struct sim *sim);
static int ext_get(FILE *f, struct sim *sim);
static int ini_get(FILE *f, struct sim *sim);
static int bnd_get(FILE *f, struct sim *sim);
