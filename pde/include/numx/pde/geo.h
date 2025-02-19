#ifndef NUMX_PDE_GEO_H
#define NUMX_PDE_GEO_H

#include <stdio.h>
#include <stdx/cap.h>
#include <stdx/cut.h>

/// @brief Geometric vertex
struct vtx {
  double x;
  double y;
  double z;

  void* ctx;
};

/// @brief Geometric line segment
struct seg {
  int vtx[2];

  void* ctx;
};

/// @brief Geometric quadrangle
struct qud {
  int vtx[4];

  void* ctx;
};

/// @brief Geometric hexahedron
struct hxd {
  int vtx[6];

  void* ctx;
};

/// @brief Complex geometric object.
struct obj {
  /// @brief Object's properties.
  struct obj_pps {
    bool with_img;
    bool with_seg;
    bool with_qud;
    bool with_hxd;
  } pps;

  /// @brief Object's vertices (required).
  struct pcut v;

  /// @brief Object's imaginary vertices.
  struct icut i;

  /// @brief Object's segments
  struct pcut s;

  /// @brief Object's quadrangles.
  struct pcut q;

  /// @brief Object's hexahedrons.
  struct pcut h;
};

int obj_new(struct obj* o);
int obj_cls(struct obj* o);

int obj_get(struct obj* o, FILE* f, struct obj_pps pps);
int obj_put(struct obj* o, FILE* f, struct obj_pps pps);

int obj_get_vtx(struct vtx* v, FILE* f);
int obj_get_seg(struct seg* s, FILE* f);
int obj_get_qud(struct qud* q, FILE* f);
int obj_get_hxd(struct hxd* h, FILE* f);

int obj_put_vtx(struct vtx* v, FILE* f);
int obj_put_seg(struct seg* s, FILE* f);
int obj_put_qud(struct qud* q, FILE* f);
int obj_put_hxd(struct hxd* h, FILE* f);

/** @brief Generate complex continuous object using supplied step function.
 *
 *  @param i initial vertex index
 *  @param s step function
 */
int obj_gen(struct obj* o, int i, struct icap s, struct obj_pps pps);

#endif  // NUMX_PDE_GEO_H