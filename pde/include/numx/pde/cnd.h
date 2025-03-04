#ifndef NUMX_PDE_CND_H
#define NUMX_PDE_CND_H

#include <numx/pde/geo.h>

typedef double (*vfun)(struct vtx*);

stdx_def_cut(ocut, vfun);

struct cnd {
  enum { DIR, NEU, ROB } type;

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

int cnd_get(struct cnd* cnd, const char* buf, struct ocut* dat);

#endif  // NUMX_PDE_CND_H