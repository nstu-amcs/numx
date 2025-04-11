#ifndef NUMX_PDE_MSH_H
#define NUMX_PDE_MSH_H

#include <numx/com/cut.h>
#include <numx/vec/vec.h>

/** Vertex. */
typedef struct vtx
{
    double x;
    double y;
    double z;
} vtx;

/** Segment. */
typedef struct seg
{
    int vtx[2];
    int pid;
} seg;

/** Quadrangle. */
typedef struct qud
{
    int vtx[4];
    int pid;
} qud;

/** Hexahedron. */
typedef struct hxd
{
    int vtx[8];
    int pid;
} hxd;

cut_def(vtx_cut, vtx);
cut_def(seg_cut, seg);
cut_def(qud_cut, qud);
cut_def(hxd_cut, hxd);

/** Unstructued mesh. */
typedef struct msh
{
    struct vtx_cut vtx; // vertices
    struct seg_cut seg; // segments
    struct qud_cut qud; // quadrangles
    struct hxd_cut hxd; // hexahedrons
} msh;

int msh_new(struct msh *msh);
int msh_cls(struct msh *msh);

/** Import mesh in Elmer format. */
int msh_imp_grd(struct msh *msh, const char *dir, const char *pfx);

/** Import mesh in CGNS format. */
int msh_exp_gns(struct msh *msh, const char *dir, const char *pfx);

/** Calculate quadrangle normal vector. */
int msh_qud_nrm(struct msh *msh, struct qud *qud, struct vec *nrm);

#endif // NUMX_PDE_MSH_H
