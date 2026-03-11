#ifndef NUMX_MSH_UMSH_H
#define NUMX_MSH_UMSH_H

#include <numx/com/cut.h>
#include <numx/vec/vec.h>

/**
 * @brief 2D Vertex.
 */
typedef struct v2d
{
    double dat[2];
} v2d;

/**
 * @brief 3D Vertex.
 */
typedef struct v3d
{
    double dat[3];
} v3d;

/**
 * @brief Generic vertex pointer.
 */
typedef union vtx_ptr
{
    struct v2d *v2d;
    struct v3d *v3d;
} vtx_ptr;

/**
 * @brief Segment.
 */
typedef struct seg
{
    int vtx[2];
    int pid; // physical group
    int qud; // encapsulating quadrangle
} seg;

/**
 * @brief Quadrangle.
 */
typedef struct qud
{
    int vtx[4];
    int pid; // physical group
    int hxd; // encapsulating hexahedron
} qud;

/**
 * @brief Hexahedron.
 */
typedef struct hxd
{
    int vtx[8];
    int pid; // physical group
} hxd;

cut_def(v2d_cut, v2d);
cut_def(v3d_cut, v3d);
cut_def(seg_cut, seg);
cut_def(qud_cut, qud);
cut_def(hxd_cut, hxd);

/**
 * @brief Unstructured mesh.
 */
typedef struct umsh
{
    enum umsh_type
    {
        MSH_C2D = 2, // cartesian 2D
        MSH_C3D = 3, // cartesian 3D
    } type;

    union
    {
        struct v2d_cut v2d;
        struct v3d_cut v3d;
    } vtx; // vertices

    struct seg_cut seg; // segments
    struct qud_cut qud; // quadrangles
    struct hxd_cut hxd; // hexahedrons
} umsh;

/**
 * @brief Segment searching function.
 *
 * For given segment, return 1 to apply the segment and 0 otherwise.
 */
typedef int (*seg_srh_fun)(struct umsh* msh, struct seg *s);

/**
 * @brief Quadrangle searching function.
 *
 * For given quadrangle, return 1 to apply the quadrangle and 0 otherwise.
 */
typedef int (*qud_srh_fun)(struct umsh* msh, struct qud *q);

int umsh_new(struct umsh *msh);
int umsh_cls(struct umsh *msh);

/**
 * @brief Import mesh from Elmer grid description.
 */
int umsh_imp_grd(struct umsh *msh, const char *dir, const char *pfx);

/**
 * @brief Import mesh from Telma grid description.
 */
int umsh_imp_tel(struct umsh *msh, const char *dir, const char *pfx);

/**
 * @brief Export mesh into CGNS.
 */
int umsh_exp_cgns(struct umsh *msh, const char *dir, const char *pfx);

/**
 * @brief Search for the segments in the mesh, creating and assigning them to the given physical
 * group.
 *
 * It's assumed that the mesh contains quadrangles as entities. For 2D mesh this is the default
 * behaviour, but for 3D mesh one must create quadrangles before searching for segments.
 */
int umsh_seg_srh(struct umsh *msh, seg_srh_fun fun, int pid);

/**
 * @brief Calculate segment's normal vector (for C2D).
 */
int umsh_seg_nrm(struct umsh *msh, struct seg *seg, struct vec *nrm);

/**
 * @brief Search for the quadrangles in the mesh, creating and assigning them to the given physical
 * group.
 *
 * It's assumed that the mesh contains hexahedrons as entities. For 3D mesh this is the default
 * behaviour.
 */
int umsh_qud_srh(struct umsh *msh, qud_srh_fun fun, int pid);

/**
 * @brief Calculate quadrangle's normal vector (for C3D).
 */
int umsh_qud_nrm(struct umsh *msh, struct qud *qud, struct vec *nrm);

/**
 * @brief Get local index in quadrangle for given global vertex.
 */
int umsh_qud_vtx_loc(struct qud *qud, int gv);

/**
 * @brief Lookup for encapsulating quadrangle for the given arbitrary vertex.
 */
int umsh_vtx_qud_lup(struct umsh *msh, union vtx_ptr vtx);

#endif // NUMX_MSH_UMSH_H
