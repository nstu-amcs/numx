#ifndef NUMX_GEO_OBJ_H
#define NUMX_GEO_OBJ_H

#include <numx/geo/geo.h>
#include <stdio.h>
#include <stdx/cut.h>

/// @brief Complex geometric object.
struct obj {
  /// @brief Object's vertices.
  struct pcut v;

  /// @brief Object's quadrangles.
  struct pcut q;

  /// @brief Object's hexahedrons.
  struct pcut h;
};

int obj_new(struct obj* o);
int obj_cls(struct obj* o);

int obj_get(struct obj* o, FILE* f);
int obj_put(struct obj* o, FILE* f);

int obj_get_vtx(struct vtx* v, FILE* f);
int obj_get_qud(struct qud* q, FILE* f);
int obj_get_hxd(struct hxd* h, FILE* f);

int obj_put_vtx(struct vtx* v, FILE* f);
int obj_put_qud(struct qud* q, FILE* f);
int obj_put_hxd(struct hxd* h, FILE* f);

/** @brief Generate complex continuous object using supplied step function.
 *
 *  @param i initial vertex index
 *  @param s step function
 */
int obj_gen_hxd(struct obj* o, int i, struct icap s);

#endif  // NUMX_GEO_OBJ_H