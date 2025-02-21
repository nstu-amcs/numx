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
    bool with_seg;
    bool with_qud;
    bool with_hxd;
  } pps;

  /// @brief Object's vertices (required).
  struct pcut v;

  /// @brief Object's segments.
  struct pcut s;

  /// @brief Object's quadrangles.
  struct pcut q;

  /// @brief Object's hexahedrons.
  struct pcut h;
};

int obj_new(struct obj* o, struct obj_pps pps);
int obj_cls(struct obj* o);

int obj_get(struct obj* o, FILE* f);
int obj_put(struct obj* o, FILE* f);

int obj_get_vtx(struct vtx* v, const char* buf);
int obj_get_seg(struct seg* s, const char* buf);
int obj_get_qud(struct qud* q, const char* buf);
int obj_get_hxd(struct hxd* h, const char* buf);

int obj_put_vtx(struct vtx* v, char* buf, int n);
int obj_put_seg(struct seg* s, char* buf, int n);
int obj_put_qud(struct qud* q, char* buf, int n);
int obj_put_hxd(struct hxd* h, char* buf, int n);

#endif  // NUMX_PDE_GEO_H