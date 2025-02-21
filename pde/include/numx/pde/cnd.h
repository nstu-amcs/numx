#ifndef NUMX_PDE_CND_H
#define NUMX_PDE_CND_H

#include <numx/geo/geo.h>

struct cnd {
  enum type { DIR, NEU, ROB } type;

  union {
    struct {
      double (*tmp)(struct vtx* v);
    } dir;

    struct {
      double (*tta)(struct vtx* v);
    } neu;

    struct {
      double (*tmp)(struct vtx* v);
      double bet;
    } rob;
  } pps;
};

#endif  // NUMX_PDE_CND_H