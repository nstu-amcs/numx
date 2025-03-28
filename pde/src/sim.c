#include <assert.h>
#include <ctype.h>
#include <dlfcn.h>
#include <stdio.h>

#include <numx/pde/sim.h>
#include <numx/pde/val.h>

cut_gen(val_cut, val, PUB);
cut_gen(mat_cut, mat, PUB);
cut_gen(obj_cut, obj, PUB);
cut_gen(bnd_cut, bnd, PUB);

cut_gen(cnd_bnd_cut, cnd_bnd, PUB);
cut_gen(cnd_ini_cut, cnd_ini, PUB);

static int get_hdr(FILE *f, struct sim *sim);
static int get_sim(FILE *f, struct sim *sim);
static int get_obj(FILE *f, struct sim *sim);
static int get_mat(FILE *f, struct sim *sim);
static int get_src(FILE *f, struct sim *sim);
static int get_ini(FILE *f, struct sim *sim);
static int get_bnd(FILE *f, struct sim *sim);

int sim_new(struct sim *sim, const char *sif)
{
    assert(sim);
    assert(sif);

    sim->ops.ell.ops.iss.ops.bcg.con.sm = 0;

    FILE *f = 0;
    int   r = 0;

    if (!(sim->msh = malloc(sizeof(struct msh)))) {
        r = -1;
        goto end;
    }

    if ((r = mat_cut_new(&sim->mat)))
        goto end;

    if ((r = val_cut_new(&sim->src)))
        goto end;

    if ((r = obj_cut_new(&sim->obj)))
        goto end;

    if ((r = bnd_cut_new(&sim->bnd)))
        goto end;

    if ((r = cnd_bnd_cut_new(&sim->cnd_bnd)))
        goto end;

    if (!(f = fopen(sif, "r"))) {
        r = -1;
        goto end;
    }

    char buf[128];

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

    if ((r = msh_new(sim->msh, sim->pps.msh.dir, sim->pps.msh.pfx)))
        goto end;
end:
    fclose(f);

    if (r) {
        free(sim->msh);

        mat_cut_cls(&sim->mat);
        val_cut_cls(&sim->src);
        obj_cut_cls(&sim->obj);
        bnd_cut_cls(&sim->bnd);

        cnd_bnd_cut_cls(&sim->cnd_bnd);
    }

    return r;
}

static int get_kv(const char *src, char *key, char *val);
static int get_str(const char **src, char *dst);
static int get_val(void *usr, const char *src, struct val *val);

static int get_hdr(FILE *f, struct sim *sim)
{
    char buf[128];
    char cmd[256];
    char str[64];

    while (fgets(buf, sizeof(buf), f)) {
        const char *cur = buf;

        while (*cur == ' ')
            ++cur;

        switch (*cur) {
            case 'E':
                return 0;
            case 'M':
                get_str(&cur, sim->pps.msh.dir);
                get_str(&cur, str);

                strcat(sim->pps.msh.dir, "/");
                strcat(sim->pps.msh.dir, str);
                strcpy(sim->pps.msh.pfx, "mesh");

                break;
            case 'I':
                get_str(&cur, str);

                if (!str[0])
                    break;

                sprintf(cmd, "gcc -o /tmp/numx_usr.so -I/usr/share/include -shared -fPIC %s", str);
                system(cmd);

                sim->pps.usr = dlopen("/tmp/numx_usr.so", RTLD_NOW);

                if (!sim->pps.usr)
                    return -1;

                break;
            case 'R':
                get_str(&cur, sim->pps.exp.dir);

                if (!sim->pps.exp.dir[0]) {
                    sim->pps.exp.dir[0] = '.';
                    sim->pps.exp.dir[1] = '\0';
                }

                break;
        }
    }

    return 0;
}

static int get_sim(FILE *f, struct sim *sim)
{
    char buf[128];
    char key[64];
    char val[64];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_kv(buf, key, val))
            return -1;

        if (!strcmp("Coordinate System", key)) {
            if (!strcmp("Cartesian 2D", val))
                sim->msh->sys = MSH_C2D;
            else if (!strcmp("Cartesian", val))
                sim->msh->sys = MSH_C3D;

            continue;
        }

        if (!strcmp("Post File", key)) {
            strcpy(strtok(val, "."), sim->pps.exp.pfx);
            strcpy(strtok(0, "."), key);

            if (!strcmp("vtu", key))
                sim->pps.exp.mod = SIM_EXP_VTU;
            else if (!strcmp("cgns", key))
                sim->pps.exp.mod = SIM_EXP_GNS;

            continue;
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
            sim->ops.ell.ops.non.err = strtod(val, 0);
            continue;
        }

        if (!strcmp("Nonlinear System Max Iterations", key)) {
            sim->ops.ell.ops.non.max = atoi(val);
            continue;
        }

        if (!strcmp("Nonlinear System Relaxation Factor", key)) {
            sim->ops.ell.ops.non.rlx = strtod(val, 0);
            continue;
        }

        if (!strcmp("Linear System Iterative Method", key)) {
            if (!strcmp("BiCGStab", val)) {
                sim->ops.ell.ops.iss.mod = ISS_BCG;
            }

            continue;
        }

        if (!strcmp("Linear System Max Iterations", key)) {
            sim->ops.ell.ops.iss.ops.bcg.ops.max = atoi(val);
            continue;
        }

        if (!strcmp("Linear System Convergence Tolerance", key)) {
            sim->ops.ell.ops.iss.ops.bcg.ops.err = strtod(val, 0);
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

        if (get_kv(buf, key, val))
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

        if (get_kv(buf, key, val))
            return -1;

        if (!strcmp("Lambda Coefficient", key)) {
            if (get_val(sim->pps.usr, val, &mat->lam))
                return -1;

            continue;
        }

        if (!strcmp("Gamma Coefficient", key)) {
            if (get_val(sim->pps.usr, val, &mat->gam))
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

        if (get_kv(buf, key, val))
            return -1;

        if (!strcmp("Field Source", key))
            if (get_val(sim->pps.usr, val, src))
                return -1;
    }

    return 0;
}

static int get_ini(FILE *f, struct sim *sim)
{
    char buf[128];
    char key[64];
    char val[64];

    if (cnd_ini_cut_dev(&sim->ops.pbc.cnd_ini, 1))
        return -1;

    struct cnd_ini *ini = &sim->ops.pbc.cnd_ini.dat[sim->ops.pbc.cnd_ini.len - 1];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_kv(buf, key, val))
            return -1;

        if (!strcmp("Field", key))
            if (get_val(sim->pps.usr, val, &ini->tgt))
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

        if (get_kv(buf, key, val))
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

            if (get_val(sim->pps.usr, val, &bnd->pps.dir.tgt))
                return -1;

            continue;
        } else if (!strcmp("Field Flux (theta)", key)) {
            bnd->type = CND_BND_NEU;

            if (get_val(sim->pps.usr, val, &bnd->pps.neu.tta))
                return -1;

        } else if (!strcmp("Robin Coefficient (beta)", key)) {
            bnd->type = CND_BND_ROB;

            if (get_val(sim->pps.usr, val, &bnd->pps.rob.bet))
                return -1;
        } else if (!strcmp("External Field", key)) {
            bnd->type = CND_BND_ROB;

            if (get_val(sim->pps.usr, val, &bnd->pps.rob.ext))
                return -1;
        }
    }

    return 0;
}

static int get_kv(const char *src, char *key, char *val)
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
