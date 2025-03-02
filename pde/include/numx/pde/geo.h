#ifndef NUMX_PDE_GEO_H
#define NUMX_PDE_GEO_H

#include <stdio.h>
#include <stdx/cap.h>
#include <stdx/cut.h>

/// @brief Geometric vertex.
typedef struct vtx {
  int n;

  double x;
  double y;
  double z;

  void* ctx;
} vtx;

/// @brief Geometric line segment.
typedef struct seg {
  int vtx[2];

  void* ctx;
} seg;

/// @brief Geometric quadrangle.
typedef struct qud {
  int vtx[4];

  void* ctx;
} qud;

/// @brief Geometric hexahedron.
typedef struct hxd {
  int vtx[8];

  void* ctx;
} hxd;

stdx_def_cut(vcut, vtx);
stdx_def_cut(scut, seg);
stdx_def_cut(qcut, qud);
stdx_def_cut(hcut, hxd);

/// @brief Complex geometric object.
struct obj {
  struct dcut ax;
  struct dcut ay;
  struct dcut az;

  /// @brief Object's vertices.
  struct vcut vtx;

  /// @brief Object's segments.
  struct scut seg;

  /// @brief Object's quadrangles.
  struct qcut qud;

  /// @brief Object's hexahedrons.
  struct hcut hxd;

  /// @brief Object's context.
  void* ctx;
};

int obj_new(struct obj* obj);
int obj_cls(struct obj* obj);

struct obj_get_ops {
  struct icap* get_seg_ctx;  // const char* buf, void** ctx
  struct icap* get_qud_ctx;  // const char* buf, void** ctx
  struct icap* get_hxd_ctx;  // const char* buf, void** ctx
};

/// @brief Get object's reference elements from the file.
int obj_get(struct obj* obj, FILE* f, struct obj_get_ops);

struct obj_put_ops {
  struct icap* put_seg_ctx;  // char* buf, size_t n, void* ctx
  struct icap* put_qud_ctx;  // char* buf, size_t n, void* ctx
  struct icap* put_hxd_ctx;  // char* buf, size_t n, void* ctx
};

/// @brief Put object's reference elements into the file.
int obj_put(struct obj* obj, FILE* f, struct obj_put_ops);

struct obj_gen_ops {
  bool with_vtx;
  bool with_seg;
  bool with_qud;
  bool with_hxd;

  struct dcap* sx;
  struct dcap* sy;
  struct dcap* sz;
};

/// @brief Generate object based on reference elements.
int obj_gen(struct obj* obj, struct obj_gen_ops ops);

#endif  // NUMX_PDE_GEO_H