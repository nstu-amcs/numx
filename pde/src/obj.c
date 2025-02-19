#include <errno.h>
#include <numx/geo/obj.h>
#include <stdio.h>
#include <stdlib.h>

int obj_new(struct obj* o) {
  if (!o) {
    errno = ENOMEM;
    return -1;
  }

  if (cut_new(&o->v))
    return -1;

  if (cut_new(&o->q))
    return -1;

  if (cut_new(&o->h))
    return -1;

  o->v.ctl = true;
  o->q.ctl = true;
  o->h.ctl = true;

  return 0;
}

int obj_cls(struct obj* o) {
  if (!o) {
    errno = ENOMEM;
    return -1;
  }

  cut_cls(&o->v);
  cut_cls(&o->q);
  cut_cls(&o->h);

  return 0;
}

int obj_get(struct obj* o, FILE* f) {
  if (!o || !f) {
    errno = EINVAL;
    return -1;
  }

  int c = 0;

  while ((c = fgetc(f)) && c != EOF)
    switch (c) {
      case (int)'v':
        struct vtx* v = malloc(sizeof(struct vtx));

        if (!v) {
          errno = ENOMEM;
          return -1;
        }

        if (obj_get_vtx(v, f)) {
          free(v);
          return -1;
        }

        if (cut_add(&o->v, v)) {
          free(v);
          return -1;
        }

        break;
      case (int)'q':
        struct qud* q = malloc(sizeof(struct qud));

        if (!q) {
          errno = ENOMEM;
          return -1;
        }

        if (obj_get_qud(q, f)) {
          free(q);
          return -1;
        }

        if (cut_add(&o->q, q)) {
          free(q);
          return -1;
        }

        break;
      case (int)'h':
        struct hxd* h = malloc(sizeof(struct hxd));

        if (!h) {
          errno = ENOMEM;
          return -1;
        }

        if (obj_get_hxd(h, f)) {
          free(h);
          return -1;
        }

        if (cut_add(&o->h, h)) {
          free(h);
          return -1;
        }

        break;
    }

  return 0;
}

int obj_get_vtx(struct vtx* v, FILE* f) {
  if (!v || !f) {
    errno = EINVAL;
    return -1;
  }

  if (fscanf(f, "%lf %lf %lf", &v->x, &v->y, &v->z) != 3)
    return -1;

  return 0;
}

int obj_get_qud(struct qud* q, FILE* f) {
  if (!q || !f) {
    errno = EINVAL;
    return -1;
  }

  if (fscanf(f, "%d %d %d %d", &q->vtx[0], &q->vtx[1], &q->vtx[2], &q->vtx[3]) != 4)
    return -1;

  return 0;
}

int obj_get_hxd(struct hxd* h, FILE* f) {
  if (!h || !f) {
    errno = EINVAL;
    return -1;
  }

  if (fscanf(f, "%d %d %d %d %d %d %d %d", &h->vtx[0], &h->vtx[1], &h->vtx[2], &h->vtx[3], &h->vtx[4], &h->vtx[5], &h->vtx[6], &h->vtx[7]) != 8)
    return -1;

  return 0;
}

int obj_put(struct obj* o, FILE* f) {
  if (!o || !f) {
    errno = EINVAL;
    return -1;
  }

  for (int i = 0; i < o->v.len; ++i) {
    if (fprintf(f, "v ") < 0)
      return -1;

    if (obj_put_vtx(o->v.dat[i], f))
      return -1;

    if (fprintf(f, "\n") < 0)
      return -1;
  }

  if (fprintf(f, "\n") < 0)
    return -1;

  for (int i = 0; i < o->q.len; ++i) {
    if (fprintf(f, "q ") < 0)
      return -1;

    if (obj_put_qud(o->q.dat[i], f))
      return -1;

    if (fprintf(f, "\n") < 0)
      return -1;
  }

  if (fprintf(f, "\n") < 0)
    return -1;

  for (int i = 0; i < o->h.len; ++i) {
    if (fprintf(f, "h ") < 0)
      return -1;

    if (obj_put_hxd(o->h.dat[i], f))
      return -1;

    if (fprintf(f, "\n") < 0)
      return -1;
  }

  return 0;
}

int obj_put_vtx(struct vtx* v, FILE* f) {
  if (!v || !f) {
    errno = EINVAL;
    return -1;
  }

  if (fprintf(f, "%lf %lf %lf", v->x, v->y, v->z) < 0)
    return -1;

  return 0;
}

int obj_put_qud(struct qud* q, FILE* f) {
  if (!q || !f) {
    errno = EINVAL;
    return -1;
  }

  if (fprintf(f, "%d %d %d %d", q->vtx[0], q->vtx[1], q->vtx[2], q->vtx[3]) < 0)
    return -1;

  return 0;
}

int obj_put_hxd(struct hxd* h, FILE* f) {
  if (!h || !f) {
    errno = EINVAL;
    return -1;
  }

  if (fprintf(f, "%d %d %d %d %d %d %d %d", h->vtx[0], h->vtx[1], h->vtx[2], h->vtx[3], h->vtx[4], h->vtx[5], h->vtx[6], h->vtx[7]) < 0)
    return -1;

  return 0;
}

int obj_gen_hxd(struct obj* o, int i, struct icap s) {
  if (!o || i < 0 || i >= o->v.len || !s.call) {
    errno = EINVAL;
    return -1;
  }

  double xs = 0;
  double ys = 0;
  double zs = 0;

  struct vtx* pv = o->v.dat[i];
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

      if (cut_add(&o->v, nv)) {
        free(nv);
        return -1;
      }

      pv = nv;

      continue;
    }

    pv = o->v.dat[rb];

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

      if (cut_add(&o->v, nv)) {
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

    pv = o->v.dat[pb];

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

      if (cut_add(&o->v, nv)) {
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