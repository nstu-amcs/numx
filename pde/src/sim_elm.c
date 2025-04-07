#include <assert.h>
#include <ctype.h>
#include <dlfcn.h>
#include <stdio.h>

#include <numx/pde/sim.h>

static int get_hdr(FILE *f, struct sim *sim);
static int get_sim(FILE *f, struct sim *sim);
static int get_obj(FILE *f, struct sim *sim);
static int get_mat(FILE *f, struct sim *sim);
static int get_src(FILE *f, struct sim *sim);
static int get_ini(FILE *f, struct sim *sim);
static int get_bnd(FILE *f, struct sim *sim);

static int get_ent(const char *src, char *key, char *val);
static int get_str(const char **src, char *dst);
static int get_val(void *usr, const char *src, struct val *val);

int sim_imp_elm(struct sim *sim, const char *sif)
{
    assert(sim);
    assert(sif);

    int   r = 0;
    FILE *f = 0;

    if (!(f = fopen(sif, "r"))) {
        r = -1;
        goto end;
    }

    char buf[256];

    while (fgets(buf, sizeof(buf), f)) {
        switch (*buf) {
            case 'H':
                if ((r = get_hdr(f, sim)))
                    goto end;

                break;
            case 'S':
                if ((r = get_sim(f, sim)))
                    goto end;

                break;
            case 'B':
                if (!strncmp(buf, "Body Force", 10)) {
                    if ((r = get_src(f, sim)))
                        goto end;

                    continue;
                }

                if (!strncmp(buf, "Boundary Condition", 18)) {
                    if ((r = get_bnd(f, sim)))
                        goto end;

                    continue;
                }

                if (!strncmp(buf, "Body", 4)) {
                    if ((r = get_obj(f, sim)))
                        goto end;

                    continue;
                }

                break;
            case 'M':
                if ((r = get_mat(f, sim)))
                    goto end;

                break;
            case 'I':
                if ((r = get_ini(f, sim)))
                    goto end;

                break;
        }
    }

end:
    fclose(f);
    return r;
}

static int get_hdr(FILE *f, struct sim *sim)
{
    char buf[256];
    char cmd[256];
    char dir[192];
    char pfx[64];
    char tmp[64];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        const char *cur = buf;

        while (*cur == ' ')
            ++cur;

        switch (*cur) {
            case 'M':
                get_str(&cur, dir);
                get_str(&cur, tmp);

                strcat(dir, "/");
                strcat(dir, tmp);
                strcpy(pfx, "mesh");

                break;
            case 'I':
                get_str(&cur, tmp);

                if (!tmp[0])
                    break;

                sprintf(cmd, "gcc -o /tmp/numx_usr.so -I/usr/share/include -shared -fPIC %s", tmp);
                system(cmd);

                sim->ops.usr = dlopen("/tmp/numx_usr.so", RTLD_NOW);

                if (!sim->ops.usr)
                    return -1;

                break;
            case 'R':
                get_str(&cur, sim->ops.exp.dir);

                if (!sim->ops.exp.dir[0]) {
                    sim->ops.exp.dir[0] = '.';
                    sim->ops.exp.dir[1] = '\0';
                }

                break;
        }
    }

    if (msh_imp_grd(sim->msh, dir, pfx))
        return -1;

    return 0;
}

static int get_sim(FILE *f, struct sim *sim)
{
    char buf[128];
    char key[64];
    char val[64];

    struct sim_ops *ops = &sim->ops;
    struct fem_ops *fem = &sim->fem->ops;

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_ent(buf, key, val))
            return -1;

        if (!strcmp("Post File", key)) {
            strcpy(ops->exp.pfx, strtok(val, "."));
            strcpy(key, strtok(0, "."));

            if (!strcmp("cgns", key)) {
                ops->exp.mod = SIM_EXP_GNS;
                ops->exp.ini.ctx = sim;
                ops->exp.ini.run = sim_exp_ini_gns;
                ops->exp.put.ctx = sim;
                ops->exp.put.run = sim_exp_put_gns;

                continue;
            }

            errno = ENOTSUP;
            return -1;
        }

        if (!strcmp("Equation", key)) {
            if (!strcmp("Elliptic Equation", val))
                sim->mod = SIM_ELL;
            else if (!strcmp("Parabolic Equation", val))
                sim->mod = SIM_PBC;
            else if (!strcmp("Hyperbolic Equation", val))
                sim->mod = SIM_HYP;

            continue;
        }

        if (!strcmp("Nonlinear System Convergence Tolerance", key)) {
            fem->non.err = strtod(val, 0);
            continue;
        }

        if (!strcmp("Nonlinear System Max Iterations", key)) {
            fem->non.max = atoi(val);

            if (fem->non.max > 1)
              fem->mod = FEM_NON;
            else
              fem->mod = FEM_STD;

            continue;
        }

        if (!strcmp("Nonlinear System Relaxation Factor", key)) {
            fem->non.rlx = strtod(val, 0);
            continue;
        }

        if (!strcmp("Linear System Iterative Method", key)) {
            if (!strcmp("BiCGStab", val)) {
                fem->iss.mod = ISS_BCG;
            }

            continue;
        }

        if (!strcmp("Linear System Max Iterations", key)) {
            fem->iss.ops.bcg.ops.max = atoi(val);
            continue;
        }

        if (!strcmp("Linear System Convergence Tolerance", key)) {
            fem->iss.ops.bcg.ops.err = strtod(val, 0);
            continue;
        }
    }

    return 0;
}

static int get_obj(FILE *f, struct sim *sim)
{
    char buf[128];
    char key[64];
    char val[64];

    if (obj_cut_dev(&sim->obj, 1))
        return -1;

    struct obj *obj = &sim->obj.dat[sim->obj.len - 1];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_ent(buf, key, val))
            return -1;

        if (!strcmp("Material", key))
            obj->mat = atoi(val) - 1;
        else if (!strcmp("Body Force", key))
            obj->src = atoi(val) - 1;
        else if (!strcmp("Initial condition", key))
            obj->ini = atoi(val) - 1;
    }

    return 0;
}

static int get_mat(FILE *f, struct sim *sim)
{
    char buf[128];
    char key[64];
    char val[64];

    if (mat_cut_dev(&sim->mat, 1))
        return -1;

    struct mat *mat = &sim->mat.dat[sim->mat.len - 1];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_ent(buf, key, val))
            return -1;

        if (!strcmp("Lambda Coefficient", key)) {
            if (get_val(sim->ops.usr, val, &mat->lam))
                return -1;

            continue;
        }

        if (!strcmp("Gamma Coefficient", key)) {
            if (get_val(sim->ops.usr, val, &mat->gam))
                return -1;

            continue;
        }
    }

    return 0;
}

static int get_src(FILE *f, struct sim *sim)
{
    char buf[128];
    char key[64];
    char val[64];

    if (val_cut_dev(&sim->src, 1))
        return -1;

    struct val *src = &sim->src.dat[sim->src.len - 1];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_ent(buf, key, val))
            return -1;

        if (!strcmp("Field Source", key))
            if (get_val(sim->ops.usr, val, src))
                return -1;
    }

    return 0;
}

static int get_ini(FILE *f, struct sim *sim)
{
    char buf[128];
    char key[64];
    char val[64];

    if (cnd_ini_cut_dev(&sim->cnd_ini, 1))
        return -1;

    struct cnd_ini *ini = &sim->cnd_ini.dat[sim->cnd_ini.len - 1];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_ent(buf, key, val))
            return -1;

        if (!strcmp("Field", key))
            if (get_val(sim->ops.usr, val, &ini->tgt))
                return -1;
    }

    return 0;
}

static int get_bnd(FILE *f, struct sim *sim)
{
    char buf[128];
    char key[64];
    char val[64];

    if (cnd_bnd_cut_dev(&sim->cnd_bnd, 1))
        return -1;

    struct cnd_bnd *bnd = &sim->cnd_bnd.dat[sim->cnd_bnd.len - 1];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_ent(buf, key, val))
            return -1;

        if (!strncmp("Target Boundaries", key, 17)) {
            for (char *ids = strtok(val, " "); ids; ids = strtok(0, " ")) {
                int idx = atoi(ids);

                if (sim->bnd.len < idx)
                    if (bnd_cut_dev(&sim->bnd, idx - sim->bnd.len))
                        return -1;

                sim->bnd.dat[idx - 1].cnd = sim->cnd_bnd.len - 1;
            }
        } else if (!strcmp("Field", key)) {
            bnd->type = CND_BND_DIR;

            if (get_val(sim->ops.usr, val, &bnd->pps.dir.tgt))
                return -1;

            continue;
        } else if (!strcmp("Field Flux (theta)", key)) {
            bnd->type = CND_BND_NEU;

            if (get_val(sim->ops.usr, val, &bnd->pps.neu.tta))
                return -1;

        } else if (!strcmp("Robin Coefficient (beta)", key)) {
            bnd->type = CND_BND_ROB;

            if (get_val(sim->ops.usr, val, &bnd->pps.rob.bet))
                return -1;
        } else if (!strcmp("External Field", key)) {
            bnd->type = CND_BND_ROB;

            if (get_val(sim->ops.usr, val, &bnd->pps.rob.ext))
                return -1;
        }
    }

    return 0;
}

static int get_ent(const char *src, char *key, char *val)
{
    while (*src == ' ')
        ++src;

    const char *beg = src;
    const char *act = src;

    while (*src != '=') {
        if (*src != ' ')
            act = src;

        ++src;
    }

    strncpy(key, beg, act - beg + 1);
    key[act - beg + 1] = 0;

    while (*(++src) == ' ') {}

    beg = src;
    act = src;

    while (*src != '\n' && *src != '\0') {
        if (*src != ' ')
            act = src;

        ++src;
    }

    strncpy(val, beg, act - beg + 1);
    val[act - beg + 1] = 0;

    return 0;
}

static int get_str(const char **src, char *dst)
{
    const char *cur = *src;
    int         len = 0;

    while (*cur != '"')
        ++cur;

    cur += 1;

    for (len = 0; *cur != '"'; ++len, ++cur)
        dst[len] = *cur;

    dst[len] = 0;
    *src = cur + 1;

    return 0;
}

static int get_val(void *usr, const char *src, struct val *val)
{
    if (isdigit(src[0])) {
        val->type = VAL_NUM;
        val->as.num = strtod(src, 0);
    } else {
        val->type = VAL_FUN;
        val->as.fun = (fun)dlsym(usr, src);

        if (!val->as.fun)
            return -1;
    }

    return 0;
}
