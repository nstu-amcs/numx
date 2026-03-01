#include <assert.h>
#include <stdio.h>

#include <numx/pde/sim.h>

int sim_imp_tel(struct sim *sim, const char *tel)
{
    assert(sim);
    assert(tel);

    int   r = 0;
    FILE *f = fopen(tel, "r");

    if (!f) {
        r = -1;
        goto end;
    }

    int n = 0;

    if (fscanf(f, "%d", &n) != 1) {
        r = -1;
        goto end;
    }

    double t = 0;
    int    pgc = 0; // physical group count
    int    pgi = 0; // physical group index

    for (int i = 0; i < n; ++i) {
        if (fscanf(f, "%lf %lf %lf %lf %lf %lf %d", &t, &t, &t, &t, &t, &t, &pgi) != 7) {
            r = -1;
            goto end;
        }

        if (pgi > pgc) {
            pgc = pgi;
        }
    }

    obj_cut_dev(&sim->obj, pgc);
    mat_cut_dev(&sim->mat, pgc);
    val_cut_dev(&sim->src, pgc);

    fseek(f, 0, SEEK_SET);
    fscanf(f, "%d", &n);

    double lam = 0;
    double src = 0;

    for (int i = 0; i < n; ++i) {
        fscanf(f, "%lf %lf %lf %lf %lf %lf %d", &t, &t, &t, &t, &lam, &src, &pgi);

        int pos = pgi - 1;

        sim->obj.dat[pos].mat = pos;
        sim->obj.dat[pos].src = pos;

        sim->mat.dat[pos].lam.type = VAL_NUM;
        sim->mat.dat[pos].lam.as.num = lam;
        sim->mat.dat[pos].gam.type = VAL_NUM;
        sim->mat.dat[pos].gam.as.num = 0;
        sim->mat.dat[pos].sig.type = VAL_NUM;
        sim->mat.dat[pos].sig.as.num = 0;
        sim->mat.dat[pos].chi.type = VAL_NUM;
        sim->mat.dat[pos].chi.as.num = 0;

        sim->src.dat[pos].type = VAL_NUM;
        sim->src.dat[pos].as.num = src;
    }

end:
    fclose(f);
    return r;
}