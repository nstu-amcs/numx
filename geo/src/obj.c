#include <errno.h>
#include <numx/geo/obj.h>
#include <stdlib.h>

int obj_gen_hxd(struct obj* o, int i, struct icap s) {
  if (!o || i < 0 || i >= o->v.len || !s.call) {
    errno = EINVAL;
    return -1;
  }

  double xs = 0;
  double ys = 0;
  double zs = 0;

  struct vtx* pv = o->v[i];
  struct vtx* nv = 0;

  int rb = i;

  double rxs = 0;
  double rys = 0;
  double rzs = 0;

  int pb = i;

  double pxs = 0;
  double pys = 0;
  double pzs = 0;

  while (1) {
    if (s.call(s.ctx, 4, pv, &xs, &ys, &zs))
      return -1;

    if (xs) {
      nv = malloc(sizeof(struct vtx));

      if (!nv) {
        errno = ENOMEM;
        return -1;
      }

      nv->x = pv->x + xs;
      nv->y = pv->y;
      nv->z = pv->z;

      if (cut_add(o->v, nv)) {
        free(nv);
        return -1;
      }

      pv = nv;

      continue;
    }

    pv = o->v[rb];

    xs = rxs;
    ys = rys;
    zs = rzs;

    if (s.call(s.ctx, 4, pv, &xs, &ys, &zs))
      return -1;

    if (ys) {
      nv = malloc(sizeof(struct vtx));

      if (!nv) {
        errno = ENOMEM;
        return -1;
      }

      nv->x = pv->x;
      nv->y = pv->y + ys;
      nv->z = pv->z;

      if (cut_add(o->v, nv)) {
        free(nv);
        return -1;
      }

      pv = nv;
      rb = o->v.len - 1;

      rxs = xs;
      rys = ys;
      rzs = zs;

      continue;
    }

    pv = o->v[pb];

    xs = pxs;
    ys = pys;
    zs = pzs;

    if (s.call(s.ctx, 4, pv, &xs, &ys, &zs))
      return -1;

    if (zs) {
      nv = malloc(sizeof(struct vtx));

      if (!nv) {
        errno = ENOMEM;
        return -1;
      }

      nv->x = pv->x;
      nv->y = pv->y;
      nv->z = pv->z + zs;

      if (cut_add(o->v, nv)) {
        free(nv);
        return -1;
      }

      pv = nv;
      pb = o->v.len - 1;

      pxs = xs;
      pys = ys;
      pzs = zs;

      continue;
    }

    break;
  }

  return 0;
}