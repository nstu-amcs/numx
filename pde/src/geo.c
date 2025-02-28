#include <errno.h>
#include <numx/pde/geo.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdx/cap.h>
#include <stdx/log.h>

stdx_gen_cut(vcut, vtx, STDX_PUB);
stdx_gen_cut(scut, seg, STDX_PUB);
stdx_gen_cut(qcut, qud, STDX_PUB);
stdx_gen_cut(hcut, hxd, STDX_PUB);

int obj_new(struct obj* obj) {
  if (!obj) {
    errno = EINVAL;
    return -1;
  }

  obj->ctx = NULL;

  if (cut_new(&obj->ax))
    return -1;

  if (cut_new(&obj->ay))
    return -1;

  if (cut_new(&obj->az))
    return -1;

  if (vcut_new(&obj->vtx))
    return -1;

  if (scut_new(&obj->seg))
    return -1;

  if (qcut_new(&obj->qud))
    return -1;

  if (hcut_new(&obj->hxd))
    return -1;

  return 0;
}

int obj_cls(struct obj* obj) {
  if (!obj) {
    errno = EINVAL;
    return -1;
  }

  cut_cls(&obj->ax);
  cut_cls(&obj->ay);
  cut_cls(&obj->az);

  vcut_cls(&obj->vtx);
  scut_cls(&obj->seg);
  qcut_cls(&obj->qud);
  hcut_cls(&obj->hxd);

  return 0;
}

static int obj_get_vtx(struct vtx* v, const char* buf) {
  if (!v || !buf) {
    errno = EINVAL;
    return -1;
  }

  int n = 0;

  if (sscanf(buf, "v %lf %lf %lf%n", &v->x, &v->y, &v->z, &n) != 3)
    return -1;

  return n;
}

static int obj_get_seg(struct seg* s, const char* buf) {
  if (!s || !buf) {
    errno = EINVAL;
    return -1;
  }

  int n = 0;

  if (sscanf(buf, "s %d %d%n", &s->vtx[0], &s->vtx[1], &n) != 2)
    return -1;

  return n;
}

static int obj_get_qud(struct qud* q, const char* buf) {
  if (!q || !buf) {
    errno = EINVAL;
    return -1;
  }

  int n = 0;

  if (sscanf(buf, "q %d %d %d %d%n", &q->vtx[0], &q->vtx[1], &q->vtx[2], &q->vtx[3], &n) != 4)
    return -1;

  return n;
}

static int obj_get_hxd(struct hxd* h, const char* buf) {
  if (!h || !buf) {
    errno = EINVAL;
    return -1;
  }

  int n = 0;

  // clang-format off
  if (sscanf(buf, "h %d %d %d %d %d %d %d %d%n", 
    &h->vtx[0], &h->vtx[1], &h->vtx[2], &h->vtx[3], 
    &h->vtx[4], &h->vtx[5], &h->vtx[6], &h->vtx[7], &n) != 8) return -1;
  // clang-format on

  return n;
}

int obj_get(struct obj* obj, FILE* f, struct obj_get_ops ops) {
  if (!obj || !f) {
    errno = EINVAL;
    return -1;
  }

  int n = 0;

  int vs = obj->vtx.len;
  int ss = obj->seg.len;
  int qs = obj->qud.len;
  int hs = obj->hxd.len;

  struct seg* s = 0;
  struct qud* q = 0;
  struct hxd* h = 0;

  char buf[128];

  while (fgets(buf, sizeof(buf), f))
    switch (buf[0]) {
      case 'v':
        if (vcut_dev(&obj->vtx, 1))
          return -1;

        if ((n = obj_get_vtx(&obj->vtx.dat[vs++], buf)) == -1)
          return -1;

        break;
      case 's':
        if (scut_dev(&obj->seg, 1))
          return -1;

        s = &obj->seg.dat[ss++];

        if ((n = obj_get_seg(s, buf)) == -1)
          return -1;

        if (ops.get_seg_ctx)
          if (ops.get_seg_ctx->call(ops.get_seg_ctx->ctx, 3, buf + n, sizeof(buf) - n, &s->ctx) == -1)
            return -1;

        break;
      case 'q':
        if (qcut_dev(&obj->qud, 1))
          return -1;

        q = &obj->qud.dat[qs++];

        if ((n = obj_get_qud(q, buf)) == -1)
          return -1;

        if (ops.get_qud_ctx)
          if (ops.get_qud_ctx->call(ops.get_qud_ctx->ctx, 3, buf + n, sizeof(buf) - n, &q->ctx) == -1)
            return -1;

        break;
      case 'h':
        if (hcut_dev(&obj->hxd, 1))
          return -1;

        h = &obj->hxd.dat[hs++];

        if ((n = obj_get_hxd(h, buf)))
          return -1;

        if (ops.get_hxd_ctx)
          if (ops.get_hxd_ctx->call(ops.get_hxd_ctx->ctx, 3, buf + n, sizeof(buf) - n, &h->ctx) == -1)
            return -1;

        break;
    }

  vcut_shr(&obj->vtx);
  scut_shr(&obj->seg);
  qcut_shr(&obj->qud);
  hcut_shr(&obj->hxd);

  return 0;
}

int obj_put(struct obj* o, FILE* f, struct obj_put_ops ops) {
  (void)ops;

  if (!o || !f) {
    errno = EINVAL;
    return -1;
  }

  errno = ENOTSUP;

  return -1;
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

static int obj_axs_div(struct dlog* a, struct dcap* s) {
  double x0 = 0;
  double x1 = 0;
  double xe = a->end->e;
  double xs = 0;

  if (log_adv(a, &x0))
    return -1;

  if (log_adv(a, &x1))
    return -1;

  while ((xs = s->call(s->ctx, 3, xs, x0, xe)) != 0) {
    double x2 = x0 + xs;

    if (x2 > x1) {
      xs = x1 - x0;
      x0 = x1;

      if (log_adv(a, &x1)) {
        if (errno = ENOENT)
          return 0;

        return -1;
      }

      continue;
    }

    x0 = x2;
    log_ins(a, x2, L);
  }

  return 0;
}

int obj_gen(struct obj* obj, struct obj_gen_ops ops) {
  if (!obj) {
    errno = EINVAL;
    return -1;
  }

  if (ops.with_seg || ops.with_qud || ops.with_hxd) {
    errno = ENOTSUP;
    return -1;
  }

  int r = 0;

  struct dlog ax;
  struct dlog ay;
  struct dlog az;

  if ((r = log_new(&ax)))
    goto end;

  if ((r = log_new(&ay)))
    goto end;

  if ((r = log_new(&az)))
    goto end;

  ax.dup = false;
  ay.dup = false;
  az.dup = false;

  ax.srt = true;
  ay.srt = true;
  az.srt = true;

  for (int i = 0; i < obj->vtx.len; ++i) {
    if (log_add(&ax, obj->vtx.dat[i].x))
      goto end;

    if (log_add(&ay, obj->vtx.dat[i].y))
      goto end;

    if (log_add(&az, obj->vtx.dat[i].z))
      goto end;
  }

  if (ops.sx && obj_axs_div(&ax, ops.sx))
    goto end;

  if (ops.sy && obj_axs_div(&ay, ops.sy))
    goto end;

  if (ops.sz && obj_axs_div(&az, ops.sz))
    goto end;

  struct vtx* vs = obj->vtx.dat;

  int i = -1;
  int j = 0;
  int l = obj->vtx.len;

  log_rst(&az);

  for (double z = 0; !log_adv(&az, &z);) {
    log_rst(&ay);

    for (double y = 0; !log_adv(&ay, &y);) {
      log_rst(&ax);

      for (double x = 0; !log_adv(&ax, &x) && j < l;) {
        i += 1;

        if (vs[j].x == x && vs[j].y == y && vs[j].z == z) {
          vs[j].n = i;
          j += 1;
        }
      }
    }
  }

  int* evx;

  for (int i = 0; i < obj->seg.len; ++i) {
    evx = obj->seg.dat[i].vtx;

    evx[0] = vs[evx[0]].n;
    evx[1] = vs[evx[1]].n;
  }

  for (int i = 0; i < obj->qud.len; ++i) {
    evx = obj->qud.dat[i].vtx;

    evx[0] = vs[evx[0]].n;
    evx[1] = vs[evx[1]].n;
    evx[2] = vs[evx[2]].n;
    evx[3] = vs[evx[3]].n;
  }

  for (int i = 0; i < obj->hxd.len; ++i) {
    evx = obj->hxd.dat[i].vtx;

    evx[0] = vs[evx[0]].n;
    evx[1] = vs[evx[1]].n;
    evx[2] = vs[evx[2]].n;
    evx[3] = vs[evx[3]].n;
    evx[4] = vs[evx[4]].n;
    evx[5] = vs[evx[5]].n;
    evx[6] = vs[evx[6]].n;
    evx[7] = vs[evx[7]].n;
  }

  if (cut_exp(&obj->ax, ax.len))
    goto end;

  if (cut_exp(&obj->ay, ay.len))
    goto end;

  if (cut_exp(&obj->az, az.len))
    goto end;

  double p = 0;

  log_rst(&ax);
  log_rst(&ay);
  log_rst(&az);

  while (!(r = log_adv(&ax, &p)))
    if ((r = cut_add(&obj->ax, p)))
      goto end;

  while (!(r = log_adv(&ay, &p)))
    if ((r = cut_add(&obj->ay, p)))
      goto end;

  while (!(r = log_adv(&az, &p)))
    if ((r = cut_add(&obj->az, p)))
      goto end;

end:
  log_cls(&ax);
  log_cls(&ay);
  log_cls(&az);

  return r;
}
