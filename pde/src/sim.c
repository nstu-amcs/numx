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
static int get_ext(FILE *f, struct sim *sim);
static int get_ini(FILE *f, struct sim *sim);
static int get_bnd(FILE *f, struct sim *sim);

int sim_new(struct sim *sim, const char *sif)
{
    assert(sim);
    assert(sif);

    sim->msh = malloc(sizeof(struct msh));

    mat_cut_new(&sim->mat);
    val_cut_new(&sim->ext);
    obj_cut_new(&sim->obj);
    bnd_cut_new(&sim->bnd);

    cnd_bnd_cut_new(&sim->cnd_bnd);

    FILE *f = fopen(sif, "r");

    if (!f)
        return -1;

    char buf[128];

    while (fgets(buf, sizeof(buf), f)) {
        switch (*buf) {
            case 'H':
                if (get_hdr(f, sim))
                    return -1;

                break;
            case 'S':
                if (get_sim(f, sim))
                    return -1;

                break;
            case 'B':
                if (!strncmp(buf, "Body Force", 10)) {
                    if (get_ext(f, sim))
                        return -1;

                    continue;
                }

                if (!strncmp(buf, "Boundary Condition", 18)) {
                    if (get_bnd(f, sim))
                        return -1;

                    continue;
                }

                if (!strncmp(buf, "Body", 4)) {
                    if (get_obj(f, sim))
                        return -1;

                    continue;
                }

                break;
            case 'M':
                if (get_mat(f, sim))
                    return -1;

                break;
            case 'I':
                if (get_ini(f, sim))
                    return -1;

                break;
        }
    }

    fclose(f);

    return 0;
}

static int get_kv(const char *src, char *key, char *val);
static int get_pt(const char **src, char *dst);

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
                get_pt(&cur, sim->pps.msh.dir);
                get_pt(&cur, str);

                strcat(sim->pps.msh.dir, "/");
                strcat(sim->pps.msh.dir, str);
                strcat(sim->pps.msh.dir, "/");

                strcpy(sim->pps.msh.pfx, "mesh");

                break;
            case 'I':
                get_pt(&cur, str);

                if (!str[0])
                    break;

                sprintf(cmd, "gcc -o /tmp/numx_usr.so -I/usr/share/include -shared -fPIC %s", str);
                system(cmd);

                sim->pps.usr = dlopen("/tmp/numx_usr.so", RTLD_NOW);

                if (!sim->pps.usr)
                    return -1;

                break;
            case 'R':
                get_pt(&cur, sim->pps.exp.dir);

                strcat(sim->pps.exp.dir, "/");

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

        get_kv(buf, key, val);

        if (!strcmp("Coordinate System", key)) {
            if (!strcmp("Cartesian 2D", val))
                sim->msh->sys = MSH_C2D;

            if (!strcmp("Cartesian", val))
                sim->msh->sys = MSH_C3D;

            continue;
        }

        if (!strcmp("Post File", key)) {
            strcpy(strtok(val, "."), sim->pps.exp.pfx);
            strcpy(strtok(0, "."), key);

            if (!strcmp("vtu", key))
                sim->pps.exp.mod = SIM_EXP_VTU;

            continue;
        }

        if (!strcmp("Equation", key)) {
            if (!strcmp("Elliptic Equation", val))
                sim->mod = SIM_ELL;

            if (!strcmp("Parabolic Equation", val))
                sim->mod = SIM_PBC;

            if (!strcmp("Hyperbolic Equation", val))
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
            if (!strcmp("BiCGStab", val))
                sim->ops.ell.ops.iss.mod = ISS_BCG;

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
    if (obj_cut_dev(&sim->obj, 1))
        return -1;

    struct obj *obj = &sim->obj.dat[sim->obj.len - 1];

    char buf[128];
    char key[64];
    char val[64];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        get_kv(buf, key, val);

        if (!strcmp("Material", key)) {
            obj->mat = atoi(val);
            continue;
        }

        if (!strcmp("Body Force", key)) {
            obj->ext = atoi(val);
            continue;
        }

        if (!strcmp("Initial condition", key)) {
            obj->ini = atoi(val);
            continue;
        }
    }

    return 0;
}

static int get_mat(FILE *f, struct sim *sim)
{
    if (mat_cut_dev(&sim->mat, 1))
        return -1;

    struct mat *mat = &sim->mat.dat[sim->mat.len - 1];

    char buf[128];
    char key[64];
    char val[64];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        get_kv(buf, key, val);

        if (!strcmp("Diffusion Coefficient", key)) {
            if (isdigit(val[0])) {
                mat->lam.type = VAL_NUM;
                mat->lam.as.num = strtod(val, 0);
            } else {
                mat->lam.type = VAL_FUN;
                mat->lam.as.fun = (fun)dlsym(sim->pps.usr, val);
            }

            continue;
        }

        if (!strcmp("Reaction Coefficient", key)) {
            if (isdigit(val[0])) {
                mat->gam.type = VAL_NUM;
                mat->gam.as.num = strtod(val, 0);
            } else {
                mat->gam.type = VAL_FUN;
                mat->gam.as.fun = (fun)dlsym(sim->pps.usr, val);
            }

            continue;
        }
    }

    return 0;
}

static int get_ext(FILE *f, struct sim *sim)
{
    if (val_cut_dev(&sim->ext, 1))
        return -1;

    struct val *ext = &sim->ext.dat[sim->ext.len - 1];

    char buf[128];
    char key[64];
    char val[64];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        get_kv(buf, key, val);

        if (!strcmp("Field Source", key)) {
            if (isdigit(val[0])) {
                ext->type = VAL_NUM;
                ext->as.num = strtod(val, 0);
            } else {
                ext->type = VAL_FUN;
                ext->as.fun = (fun)dlsym(sim->pps.usr, val);
            }

            continue;
        }
    }

    return 0;
}

static int get_ini(FILE *f, struct sim *sim)
{
    if (cnd_ini_cut_dev(&sim->ops.pbc.cnd_ini, 1))
        return -1;

    struct cnd_ini *ini = &sim->ops.pbc.cnd_ini.dat[sim->ops.pbc.cnd_ini.len - 1];

    char buf[128];
    char key[64];
    char val[64];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        get_kv(buf, key, val);

        if (!strcmp("Field", key)) {
            if (isdigit(val[0])) {
                ini->tgt.type = VAL_NUM;
                ini->tgt.as.num = strtod(val, 0);
            } else {
                ini->tgt.type = VAL_FUN;
                ini->tgt.as.fun = (fun)dlsym(sim->pps.usr, val);
            }

            continue;
        }
    }

    return 0;
}

static int get_bnd(FILE *f, struct sim *sim)
{
    if (cnd_bnd_cut_dev(&sim->cnd_bnd, 1))
        return -1;

    struct cnd_bnd *bnd = &sim->cnd_bnd.dat[sim->cnd_bnd.len - 1];

    char buf[128];
    char key[64];
    char val[64];

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        get_kv(buf, key, val);

        if (!strncmp("Target Boundaries", key, 17)) {
            char *ids = strtok(val, " ");

            while (ids) {
                int idx = atoi(ids);

                if (sim->bnd.len < idx)
                    bnd_cut_dev(&sim->bnd, idx - sim->bnd.len);

                sim->bnd.dat[idx - 1].cnd = sim->cnd_bnd.len - 1;

                ids = strtok(0, " ");
            }

            continue;
        }

        if (!strcmp("Field", key)) {
            bnd->type = CND_BND_DIR;

            if (isdigit(val[0])) {
                bnd->pps.dir.tgt.type = VAL_NUM;
                bnd->pps.dir.tgt.as.num = strtod(val, 0);
            } else {
                bnd->pps.dir.tgt.type = VAL_FUN;
                bnd->pps.dir.tgt.as.fun = (fun)dlsym(sim->pps.usr, val);
            }

            continue;
        }

        if (!strcmp("Field Flux", key)) {
            bnd->type = CND_BND_NEU;

            if (isdigit(val[0])) {
                bnd->pps.neu.tta.type = VAL_NUM;
                bnd->pps.neu.tta.as.num = strtod(val, 0);
            } else {
                bnd->pps.neu.tta.type = VAL_FUN;
                bnd->pps.neu.tta.as.fun = (fun)dlsym(sim->pps.usr, val);
            }

            continue;
        }

        if (!strcmp("Robin Coefficient", key)) {
            bnd->type = CND_BND_ROB;

            if (isdigit(val[0])) {
                bnd->pps.rob.bet.type = VAL_NUM;
                bnd->pps.rob.bet.as.num = strtod(val, 0);
            } else {
                bnd->pps.rob.bet.type = VAL_FUN;
                bnd->pps.rob.bet.as.fun = (fun)dlsym(sim->pps.usr, val);
            }

            continue;
        }

        if (!strcmp("External Field", key)) {
            bnd->type = CND_BND_ROB;

            if (isdigit(val[0])) {
                bnd->pps.rob.src.type = VAL_NUM;
                bnd->pps.rob.src.as.num = strtod(val, 0);
            } else {
                bnd->pps.rob.src.type = VAL_FUN;
                bnd->pps.rob.src.as.fun = (fun)dlsym(sim->pps.usr, val);
            }

            continue;
        }
    }

    return 0;
}

static int get_pt(const char **src, char *dst)
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
