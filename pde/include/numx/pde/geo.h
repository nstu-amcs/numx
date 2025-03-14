#ifndef NUMX_PDE_GEO_H
#define NUMX_PDE_GEO_H

#include <stdio.h>
#include <stdx/cap.h>
#include <stdx/cut.h>

/** Normal direction (should be replaced with vector-based normal calculation). */
enum norm
{
    NORM_U,
    NORM_D,
    NORM_L,
    NORM_R,
    NORM_F,
    NORM_B
};

/** Geometric vertex. */
typedef struct vtx
{
    double x;
    double y;
    double z;

    int   n;
    void *ctx;
} vtx;

/** Geometric line segment. */
typedef struct seg
{
    int vtx[2];

    void *ctx;
} seg;

/** Geometric quadrangle. */
typedef struct qud
{
    int vtx[4];

    void *ctx;
} qud;

/** Geometric hexahedron. */
typedef struct hxd
{
    int vtx[8];

    void *ctx;
} hxd;

stdx_def_cut(vtx_cut, vtx);
stdx_def_cut(seg_cut, seg);
stdx_def_cut(qud_cut, qud);
stdx_def_cut(hxd_cut, hxd);

/** Complex geometric object. */
struct obj
{
    struct dcut ax; // X-axis
    struct dcut ay; // Y-axis
    struct dcut az; // Z-axis

    /** Reference vertices. */
    struct vtx_cut vtx;

    /** Reference segments. */
    struct seg_cut seg;

    /** Reference quadrangles. */
    struct qud_cut qud;

    /** Reference hexahedrons. */
    struct hxd_cut hxd;

    void *ctx;
};

int obj_new(struct obj *obj);
int obj_cls(struct obj *obj);

struct obj_get_ops
{
    struct icap *get_vtx_ctx;
    struct icap *get_seg_ctx;
    struct icap *get_qud_ctx;
    struct icap *get_hxd_ctx;
};

/** Get object's reference elements from the file. */
int obj_get(struct obj *obj, FILE *f, struct obj_get_ops);

struct obj_gen_ops
{
    struct dcap *sx; // X-step
    struct dcap *sy; // Y-step
    struct dcap *sz; // Z-step

    double eps; // Reference point restricted neighborhood
};

/** Generate object based on reference elements. */
int obj_gen(struct obj *obj, struct obj_gen_ops ops);

/** Lookup global number for given point. */
int obj_srh(struct obj *obj, double x, double y, double z);

/** Discover normal direction for given quadrangle. */
enum norm qud_norm(struct obj *obj, struct qud *qud);

#endif // NUMX_PDE_GEO_H
