#include <assert.h>
#include <ctype.h>
#include <dlfcn.h>
#include <stdio.h>

#include <numx/pde/sim.h>

static int get_sim(FILE *f, struct sim *sim);
static int get_obj(FILE *f, struct sim *sim);
static int get_mat(FILE *f, struct sim *sim);
static int get_src(FILE *f, struct sim *sim);
static int get_ini(FILE *f, struct sim *sim);
static int get_bnd(FILE *f, struct sim *sim);

static int get_ent(const char *src, char *key, char *val);
static int get_val(struct sim *sim, char *src, struct val *val);

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

    mat_cut_shr(&sim->mat);
    val_cut_shr(&sim->src);
    obj_cut_shr(&sim->obj);
    bnd_cut_shr(&sim->bnd);
    cnd_ini_cut_shr(&sim->cnd_ini);
    cnd_bnd_cut_shr(&sim->cnd_bnd);

end:
    if (r)
        free(sim->slv);

    fclose(f);

    return r;
}

static int get_sim(FILE *f, struct sim *sim)
{
    char buf[128];
    char key[64];
    char val[64];

    struct fem_ops *fem = &((struct fem *)sim->slv)->ops;

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_ent(buf, key, val))
            return -1;

        if (!strncmp("Timestep intervals", key, 18)) {
            sim->ops.tdd.num = atoi(val);
            continue;
        }

        if (!strncmp("Timestep Sizes", key, 14)) {
            sscanf(val, "%lf %lf", &sim->ops.tdd.beg, &sim->ops.tdd.hop);
            continue;
        }

        if (!strcmp("BDF Order", key)) {
            sim->slv->ops.tdd = atoi(val);
            continue;
        }

        if (!strcmp("Equation", key)) {
            if (!strcmp("Elliptic Equation", val))
                sim->eqn = SIM_ELL;
            else if (!strcmp("Parabolic Equation", val))
                sim->eqn = SIM_PBC;
            else if (!strcmp("Hyperbolic Equation", val))
                sim->eqn = SIM_HYP;

            continue;
        }

        if (!strcmp("Solution Mode", key)) {
            if (!strcmp("Standard", val))
                sim->mod = SIM_STD;
            else if (!strcmp("Harmonic", val))
                sim->mod = SIM_HMC;

            continue;
        }

        if (!strcmp("Basis", key)) {
            if (!strcmp("Linear", val))
                fem->bss = FEM_BSS_LIN;

            continue;
        }

        if (!strcmp("Nonlinear System Convergence Tolerance", key)) {
            sim->slv->ops.non.ops.err = strtod(val, 0);
            continue;
        }

        if (!strcmp("Nonlinear System Max Iterations", key)) {
            sim->slv->ops.non.ops.max = atoi(val);
            continue;
        }

        if (!strcmp("Linear System Iterative Method", key)) {
            if (!strcmp("BiCGStab", val)) {
                sim->slv->ops.iss.mod = ISS_BCG;
            }

            continue;
        }

        if (!strcmp("Linear System Max Iterations", key)) {
            sim->slv->ops.iss.ops.bcg.ops.max = atoi(val);
            continue;
        }

        if (!strcmp("Linear System Convergence Tolerance", key)) {
            sim->slv->ops.iss.ops.bcg.ops.err = strtod(val, 0);
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

    mat->lam.type = VAL_NUM;
    mat->lam.as.num = 0;
    mat->gam.type = VAL_NUM;
    mat->gam.as.num = 0;
    mat->sig.type = VAL_NUM;
    mat->sig.as.num = 0;
    mat->chi.type = VAL_NUM;
    mat->chi.as.num = 0;

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_ent(buf, key, val))
            return -1;

        if (!strcmp("Lambda Coefficient", key)) {
            if (get_val(sim, val, &mat->lam))
                return -1;

            if (mat->lam.ops.dep)
                NON_LAM_SET(sim->slv->ops.non.map);

            continue;
        }

        if (!strcmp("Gamma Coefficient", key)) {
            if (get_val(sim, val, &mat->gam))
                return -1;

            if (mat->gam.ops.dep)
                NON_GAM_SET(sim->slv->ops.non.map);

            continue;
        }

        if (!strcmp("Sigma Coefficient", key)) {
            if (get_val(sim, val, &mat->sig))
                return -1;

            if (mat->sig.ops.dep)
                NON_SIG_SET(sim->slv->ops.non.map);

            continue;
        }

        if (!strcmp("Chi Coefficient", key)) {
            if (get_val(sim, val, &mat->chi))
                return -1;

            if (mat->chi.ops.dep)
                NON_CHI_SET(sim->slv->ops.non.map);

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

    src->type = VAL_NUM;
    src->as.num = 0;

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_ent(buf, key, val))
            return -1;

        if (!strcmp("Field Source", key)) {
            if (get_val(sim, val, src))
                return -1;

            if (src->ops.dep)
                NON_SRC_SET(sim->slv->ops.non.map);
        }
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

    ini->tgt.type = VAL_NUM;
    ini->tgt.as.num = 0;

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == 'E')
            break;

        if (get_ent(buf, key, val))
            return -1;

        if (!strcmp("Field", key))
            if (get_val(sim, val, &ini->tgt))
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

            if (get_val(sim, val, &bnd->pps.dir.tgt))
                return -1;

            continue;
        } else if (!strcmp("Field Flux (theta)", key)) {
            bnd->type = CND_BND_NEU;

            if (get_val(sim, val, &bnd->pps.neu.tta))
                return -1;

            if (bnd->pps.neu.tta.ops.dep)
                NON_TTA_SET(sim->slv->ops.non.map);

            continue;
        } else if (!strcmp("Robin Coefficient (beta)", key)) {
            bnd->type = CND_BND_ROB;

            if (get_val(sim, val, &bnd->pps.rob.bet))
                return -1;

            if (bnd->pps.rob.bet.ops.dep)
                NON_BET_SET(sim->slv->ops.non.map);

            continue;
        } else if (!strcmp("External Field", key)) {
            bnd->type = CND_BND_ROB;

            if (get_val(sim, val, &bnd->pps.rob.ext))
                return -1;

            if (bnd->pps.rob.ext.ops.dep)
                NON_EXT_SET(sim->slv->ops.non.map);

            continue;
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

static int get_val(struct sim *sim, char *src, struct val *val)
{
    char *fun = strtok(src, ";");

    if (isdigit(fun[0])) {
        val->type = VAL_NUM;
        val->as.num = strtod(fun, 0);
    } else {
        val->type = VAL_FUN;
        val->as.fun = dlsym(sim->ops.usr.hdl, fun);

        if (!val->as.fun)
            return -1;
    }

    val->ops.dep = false;
    val->ops.dif = NULL;

    if ((fun = strtok(0, ";"))) {
        val->ops.dep = true;

        if (strcmp("num", fun))
            val->ops.dif = dlsym(sim->ops.usr.hdl, fun);
    }

    return 0;
}
