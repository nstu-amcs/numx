#include <assert.h>
#include <stdio.h>

#include <numx/pde/msh.h>

static int get_vtx(struct msh *msh, FILE *f);

static int get_hxd_ems(struct msh *msh, FILE *f);
static int get_qud_bnd(struct msh *msh, FILE *f);

int msh_imp_grd(struct msh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);

    int r = 0;
    int n = 0;
    int e = 0;
    int b = 0;

    FILE *hdr = 0;
    FILE *vtx = 0;
    FILE *ems = 0;
    FILE *bnd = 0;

    char path[64];

    sprintf(path, "%s/%s.header", dir, pfx);

    if (!(hdr = fopen(path, "r"))) {
        r = -1;
        goto end;
    }

    sprintf(path, "%s/%s.nodes", dir, pfx);

    if (!(vtx = fopen(path, "r"))) {
        r = -1;
        goto end;
    }

    sprintf(path, "%s/%s.elements", dir, pfx);

    if (!(ems = fopen(path, "r"))) {
        r = -1;
        goto end;
    }

    sprintf(path, "%s/%s.boundary", dir, pfx);

    if (!(bnd = fopen(path, "r"))) {
        r = -1;
        goto end;
    }

    if (fscanf(hdr, "%d %d %d", &n, &e, &b) != 3) {
        r = -1;
        goto end;
    }
    if ((r = vec_cut_dev(&msh->vtx, n)))
        goto end;

    if ((r = hxd_cut_dev(&msh->hxd, e)))
        goto end;

    if ((r = qud_cut_dev(&msh->qud, b)))
        goto end;

    if ((r = get_vtx(msh, vtx)))
        goto end;

    if ((r = get_qud_bnd(msh, bnd)))
        goto end;

    if ((r = get_hxd_ems(msh, ems)))
        goto end;

    vec_cut_shr(&msh->vtx);
    seg_cut_shr(&msh->seg);
    qud_cut_shr(&msh->qud);
    hxd_cut_shr(&msh->hxd);

end:
    fclose(hdr);
    fclose(vtx);
    fclose(ems);
    fclose(bnd);

    return r;
}

int msh_exp_grd(struct msh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);

    errno = ENOTSUP;

    return -1;
}

// clang-format off

static int get_vtx(struct msh *msh, FILE *f)
{

    for (int i = 0, j; i < msh->vtx.len; ++i) {
        struct vec *vtx = &msh->vtx.dat[i];

        if (vec_new(vtx, 3))
            return -1;

        if (fscanf(f, "%d %d %lf %lf %lf", &j, &j, 
              &vtx[i].dat[0], 
              &vtx[i].dat[1], 
              &vtx[i].dat[2]) != 5)
            return -1;
    }

    return 0;
}

static int get_hxd_ems(struct msh *msh, FILE *f)
{
    struct hxd *hxd = msh->hxd.dat;

    for (int i = 0, j; i < msh->hxd.len; ++i) {
        if (fscanf(f, "%d %d %d %d %d %d %d %d %d %d %d", &j, &hxd[i].pid, &j, 
              &hxd[i].vtx[4], 
              &hxd[i].vtx[0],
              &hxd[i].vtx[2],
              &hxd[i].vtx[6], 
              &hxd[i].vtx[5], 
              &hxd[i].vtx[1], 
              &hxd[i].vtx[3], 
              &hxd[i].vtx[7]) != 11)
            return -1;

        hxd[i].pid -= 1;
        hxd[i].vtx[0] -= 1;
        hxd[i].vtx[1] -= 1;
        hxd[i].vtx[2] -= 1;
        hxd[i].vtx[3] -= 1;
        hxd[i].vtx[4] -= 1;
        hxd[i].vtx[5] -= 1;
        hxd[i].vtx[6] -= 1;
        hxd[i].vtx[7] -= 1;
    }

    return 0;
}

static int qud_srt(struct msh* msh, struct qud* qud);

static int get_qud_bnd(struct msh *msh, FILE *f)
{
    struct qud *qud = msh->qud.dat;

    for (int i = 0, j; i < msh->qud.len; ++i) {
        if (fscanf(f, "%d %d %d %d %d %d %d %d %d", &j, &qud[i].pid, &qud[i].hxd, &j, &j, 
              &qud[i].vtx[0], 
              &qud[i].vtx[1],
              &qud[i].vtx[2], 
              &qud[i].vtx[3]) != 9)
            return -1;

        qud[i].pid -= 1;
        qud[i].hxd -= 1;
        qud[i].vtx[0] -= 1;
        qud[i].vtx[1] -= 1;
        qud[i].vtx[2] -= 1;
        qud[i].vtx[3] -= 1;

        qud_srt(msh, &qud[i]);
    }

    return 0;
}

static inline int cmp(struct msh* msh, int a, int b)
{
    double av = msh->vtx.dat[a].dat[2];
    double bv = msh->vtx.dat[b].dat[2];

    if (av < bv)
      return 0;

    if (bv < av)
      return 1;

    av = msh->vtx.dat[a].dat[1];
    bv = msh->vtx.dat[b].dat[1];

    if (av < bv)
      return 0;

    if (bv < av)
      return 1;

    av = msh->vtx.dat[a].dat[0];
    bv = msh->vtx.dat[b].dat[0];

    if (av < bv)
      return 0;

    return 1;
}

static int qud_srt(struct msh* msh, struct qud* qud)
{
    int t;
    int *v = qud->vtx;

    if (cmp(msh, v[0], v[1])) {
      t = v[0];
      v[0] = v[1];
      v[1] = t;
    }

    if (cmp(msh, v[2], v[3])) {
      t = v[2];
      v[2] = v[3];
      v[3] = t;
    }

    if (cmp(msh, v[0], v[2])) {
      t = v[0];
      v[0] = v[2];
      v[2] = t;

      if (cmp(msh, v[2], v[3])) {
        t = v[1];
        v[1] = v[3];
        v[3] = t;
      }
    } else {
      if (cmp(msh, v[1], v[2])) {
        t = v[1];
        v[1] = v[2];
        v[2] = t;

        if (cmp(msh, v[2], v[3])) {
          t = v[2];
          v[2] = v[3];
          v[3] = t;
        }
      }
    }

    return 0;
}
