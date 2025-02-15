#ifndef NUMX_GEO_H
#define NUMX_GEO_H

struct vtx {
  double x;
  double y;
  double z;
};

struct rct {
  int vtx[4];
};

struct hxd {
  int vtx[6];
};

#endif  // NUMX_GEO_H