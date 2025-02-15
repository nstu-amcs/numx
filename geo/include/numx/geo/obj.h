#ifndef NUMX_GEO_OBJ_H
#define NUMX_GEO_OBJ_H

#include <numx/geo/geo.h>
#include <stdio.h>

struct obj {
  struct vtx** v;
  struct rct** r;
  struct hxd** h;
};

int obj_get(struct obj* o, FILE* f);
int obj_put(struct obj* o, FILE* f);

#endif  // NUMX_GEO_OBJ_H