#ifndef NUMX_PDE_MSH_H
#define NUMX_PDE_MSH_H

#include <numx/com/cut.h>

typedef struct vtx
{
    double x;
    double y;
    double z;
} vtx;

typedef struct seg
{
    int vtx[2];
    int pid;
} seg;

typedef struct qud
{
    int vtx[4];
    int pid;
} qud;

typedef struct hxd
{
    int vtx[8];
    int pid;
} hxd;

cut_def(vtx_cut, vtx);
cut_def(seg_cut, seg);
cut_def(qud_cut, qud);
cut_def(hxd_cut, hxd);

typedef struct msh
{
    enum
    {
        MSH_C2D,
        MSH_C3D,
    } sys;

    struct vtx_cut vtx;
    struct seg_cut seg;
    struct qud_cut qud;
    struct hxd_cut hxd;
} msh;

int msh_new(struct msh *msh, const char *dir);
int msh_cls(struct msh *msh);

#endif // NUMX_PDE_MSH_H
