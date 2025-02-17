#ifndef NUMX_GEO_GEO_H
#define NUMX_GEO_GEO_H

/// @brief Geometric vertex
struct vtx {
  double x;
  double y;
  double z;

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

#endif  // NUMX_GEO_GEO_H