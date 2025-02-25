#include <errno.h>
#include <numx/pde/geo.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdx/log.h>

int obj_new(struct obj* o, struct obj_pps p) {
  if (!o) {
    errno = ENOMEM;
    return -1;
  }

  o->pps = p;

  if (cut_new(&o->v))
    return -1;

  if (cut_new(&o->s))
    return -1;

  if (cut_new(&o->q))
    return -1;

  if (cut_new(&o->h))
    return -1;

  o->v.ctl = true;
  o->s.ctl = true;
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
  cut_cls(&o->s);
  cut_cls(&o->q);
  cut_cls(&o->h);

  return 0;
}

int obj_get(struct obj* o, FILE* f) {
  if (!o || !f) {
    errno = EINVAL;
    return -1;
  }

  char buf[64];

  struct vtx* v = 0;
  struct seg* s = 0;
  struct qud* q = 0;
  struct hxd* h = 0;

  while (fgets(buf, sizeof(buf), f))
    switch (buf[0]) {
      case 'v':
        v = malloc(sizeof(struct vtx));

        if (!v) {
          errno = ENOMEM;
          return -1;
        }

        if (obj_get_vtx(v, buf)) {
          free(v);
          return -1;
        }

        if (cut_add(&o->v, v)) {
          free(v);
          return -1;
        }

        break;
      case 's':
        if (!o->pps.with_seg)
          continue;

        s = malloc(sizeof(struct seg));

        if (!s) {
          errno = ENOMEM;
          return -1;
        }

        if (obj_get_seg(s, buf)) {
          free(s);
          return -1;
        }

        if (cut_add(&o->s, s)) {
          free(s);
          return -1;
        }

        break;
      case 'q':
        if (!o->pps.with_qud)
          continue;

        q = malloc(sizeof(struct qud));

        if (!q) {
          errno = ENOMEM;
          return -1;
        }

        if (obj_get_qud(q, buf)) {
          free(q);
          return -1;
        }

        if (cut_add(&o->q, q)) {
          free(q);
          return -1;
        }

        break;
      case 'h':
        if (!o->pps.with_hxd)
          continue;

        h = malloc(sizeof(struct hxd));

        if (!h) {
          errno = ENOMEM;
          return -1;
        }

        if (obj_get_hxd(h, buf)) {
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

int obj_get_vtx(struct vtx* v, const char* buf) {
  if (!v || !buf) {
    errno = EINVAL;
    return -1;
  }

  if (sscanf(buf, "v %lf %lf %lf", &v->x, &v->y, &v->z) != 3)
    return -1;

  return 0;
}

int obj_get_seg(struct seg* s, const char* buf) {
  if (!s || !buf) {
    errno = EINVAL;
    return -1;
  }

  if (sscanf(buf, "s %d %d", &s->vtx[0], &s->vtx[1]) != 2)
    return -1;

  return 0;
}

int obj_get_qud(struct qud* q, const char* buf) {
  if (!q || !buf) {
    errno = EINVAL;
    return -1;
  }

  if (sscanf(buf, "q %d %d %d %d", &q->vtx[0], &q->vtx[1], &q->vtx[2], &q->vtx[3]) != 4)
    return -1;

  return 0;
}

int obj_get_hxd(struct hxd* h, const char* buf) {
  if (!h || !buf) {
    errno = EINVAL;
    return -1;
  }

  if (sscanf(buf, "h %d %d %d %d %d %d %d %d", &h->vtx[0], &h->vtx[1], &h->vtx[2], &h->vtx[3], &h->vtx[4], &h->vtx[5], &h->vtx[6], &h->vtx[7]) != 8)
    return -1;

  return 0;
}

int obj_put(struct obj* o, FILE* f) {
  if (!o || !f) {
    errno = EINVAL;
    return -1;
  }

  char buf[64];

  for (int i = 0; i < o->v.len; ++i) {
    if (obj_put_vtx(o->v.dat[i], buf, sizeof(buf)))
      return -1;

    if (fputs(buf, f) == EOF)
      return -1;

    if (fputc('\n', f) == EOF)
      return -1;
  }

  if (o->pps.with_seg) {
    if (fputc('\n', f) == EOF)
      return -1;

    for (int i = 0; i < o->s.len; ++i) {
      if (obj_put_seg(o->s.dat[i], buf, sizeof(buf)))
        return -1;

      if (fputs(buf, f) == EOF)
        return -1;

      if (fputc('\n', f) == EOF)
        return -1;
    }
  }

  if (o->pps.with_qud) {
    if (fputc('\n', f) == EOF)
      return -1;

    for (int i = 0; i < o->q.len; ++i) {
      if (obj_put_qud(o->q.dat[i], buf, sizeof(buf)))
        return -1;

      if (fputs(buf, f) == EOF)
        return -1;

      if (fputc('\n', f) == EOF)
        return -1;
    }
  }

  if (o->pps.with_hxd) {
    if (fputc('\n', f) == EOF)
      return -1;

    for (int i = 0; i < o->h.len; ++i) {
      if (obj_put_seg(o->h.dat[i], buf, sizeof(buf)))
        return -1;

      if (fputs(buf, f) == EOF)
        return -1;

      if (fputc('\n', f) == EOF)
        return -1;
    }
  }

  return 0;
}

int obj_put_vtx(struct vtx* v, char* buf, int n) {
  if (!v || !buf) {
    errno = EINVAL;
    return -1;
  }

  int r = snprintf(buf, n, "v %lf %lf %lf", v->x, v->y, v->z);

  if (r < 0 || r >= n)
    return -1;

  return 0;
}

int obj_put_seg(struct seg* s, char* buf, int n) {
  if (!s || !buf) {
    errno = EINVAL;
    return -1;
  }

  int r = snprintf(buf, n, "s %d %d", s->vtx[0], s->vtx[1]);

  if (r < 0 || r >= n)
    return -1;

  return 0;
}

int obj_put_qud(struct qud* q, char* buf, int n) {
  if (!q || !buf) {
    errno = EINVAL;
    return -1;
  }

  int r = snprintf(buf, n, "q %d %d %d %d", q->vtx[0], q->vtx[1], q->vtx[2], q->vtx[3]);

  if (r < 0 || r >= n)
    return -1;

  return 0;
}

int obj_put_hxd(struct hxd* h, char* buf, int n) {
  if (!h || !buf) {
    errno = EINVAL;
    return -1;
  }

  int r = snprintf(buf, n, "h %d %d %d %d %d %d %d %d", h->vtx[0], h->vtx[1], h->vtx[2], h->vtx[3], h->vtx[4], h->vtx[5], h->vtx[6], h->vtx[7]);

  if (r < 0 || r >= n)
    return -1;

  return 0;
}

// clang-format off
int obj_gen_img(struct obj* o, 
  struct dcap* xs,
  struct dcap* ys,
  struct dcap* zs,
  struct vcap* m
) {
  // clang-format on

  int r = 0;

  if (!o) {
    errno = EINVAL;
    return -1;
  }

  struct vtx** vv = (struct vtx**)o->v.dat;

  struct plog ls;
  struct ilog xa;
  struct ilog ya;

  if ((r = log_new(&ls)))
    goto end;

  if ((r = log_new(&xa)))
    goto end;

  if ((r = log_new(&ya)))
    goto end;

  xa.srt = true;
  ya.srt = true;

  struct plog* l = 0;
  double y = vv[0]->y - 1;

  for (int i = 0; i < o->v.len; ++i) {
    struct vtx* v = vv[i];

    if (v->y != y) {
      l = malloc(sizeof(struct plog));

      if (!l) {
        errno = ENOMEM;
        r = -1;
        goto end;
      }

      if ((r = log_new(l))) {
        free(l);
        goto end;
      }

      if ((r = log_add(&ls, l))) {
        log_cls(l);
        free(l);
        goto end;
      }
    }

    if ((r = log_add(l, v)))
      goto end;

    if ((r = log_add(&xa, v->x)) && errno != EALREADY)
      goto end;

    if ((r = log_add(&ya, v->y)) && errno != EALREADY)
      goto end;
  }

  if (xs) {
    double x = 0;
    double s = 0;
    double n = 0;
    double p = 0;

    log_adv(&xa, &x);

    while (log_adv(&xa, &n) == 1) {
      s = xs->call(xs->ctx, 2, x, s);
      n = x + s;

      
    }
  }

end:
  for (struct plog* i = 0; log_adv(&ls, &i);) {
    log_cls(i);
    free(i);
  }

  log_cls(&ls);
  log_cls(&xa);
  log_cls(&ya);

  return r;
}
