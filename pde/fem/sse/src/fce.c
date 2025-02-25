#include "const.h"

#include <fem/sse/fce.h>
#include <vec/mtx.h>

#include <errno.h>
#include <stdlib.h>

int fce_ini(struct fce** h) {
  if (!h) {
    errno = EINVAL;
    return -1;
  }

  struct fce* f = malloc(sizeof(struct fce));

  if (!f) {
    errno = ENOMEM;
    return -1;
  }

  f->vtx[0] = 0;
  f->vtx[1] = 0;
  f->vtx[2] = 0;
  f->vtx[3] = 0;

  f->cnd.pps.rob.bet = 0;
  f->cnd.pps.rob.tmp = 0;

  f->loc.rob.b = 0;
  f->loc.rob.m = 0;

  *h = f;

  return 0;
}

int fce_cls(struct fce** h) {
  if (!h || !(*h)) {
    errno = EINVAL;
    return -1;
  }

  struct fce* f = *h;

  switch (f->cnd.type) {
    case DIR:
      break;
    case ROB:
      if (f->loc.rob.b)
        vec_cls(f->loc.rob.b);

      if (f->loc.rob.m)
        mtx_cls(f->loc.rob.m);

      break;
    case NEU:
      if (f->loc.neu.b)
        vec_cls(f->loc.neu.b);

      break;
  }

  free(f);
  *h = 0;

  return 0;
}

int fce_new(struct fce* f) {
  if (!f) {
    errno = EINVAL;
    return -1;
  }

  switch (f->cnd.type) {
    case DIR:
      return 0;
    case ROB:
      f->loc.rob.b = vec_new(4);
      f->loc.rob.m = mtx_new(4);
      return 0;
    case NEU:
      f->loc.neu.b = vec_new(4);
      return 0;
  }

  return 0;
}

int fce_sget(struct fce* f, const char* buf, double (**fun)(struct vtx*)) {
  if (buf[0] != 'f') {
    errno = EINVAL;
    return -1;
  }

  int pos, pty;

  buf += 1;

  sscanf(buf, "%d |%n", &pty, &pos);
  buf += pos;

  f->cnd.type = pty - 1;

  for (int j = 0; j < 4; ++j) {
    sscanf(buf, "%d%n", &pty, &pos);
    buf += pos;

    f->vtx[j] = pty - 1;
  }

  switch (f->cnd.type) {
    case DIR:
      sscanf(buf, " | %d%n", &pty, &pos);
      buf += pos;

      f->cnd.pps.dir.tmp = fun[pty - 1];

      break;
    case NEU:
      sscanf(buf, " | %d%n", &pty, &pos);
      buf += pos;

      f->cnd.pps.neu.tta = fun[pty - 1];

      break;
    case ROB:
      double bet;

      sscanf(buf, " | %d | %lf%n", &pty, &bet, &pos);
      buf += pos;

      f->cnd.pps.rob.tmp = fun[pty - 1];
      f->cnd.pps.rob.bet = bet;

      break;
  }

  fce_new(f);

  return 0;
}

static int fce_nrm(struct fce* f, struct vtx** v) {
  struct vtx* a = v[f->vtx[0]];
  struct vtx* b = v[f->vtx[1]];
  struct vtx* c = v[f->vtx[2]];

  if (a->x == c->x && b->x == c->x)
    return 0;

  if (a->y == c->y && b->y == c->y)
    return 1;

  return 2;
}

int fce_evo(struct fce* f, struct vtx** v) {
  if (f->cnd.type == DIR)
    return 0;

  double mxi[2][2];
  double mzt[2][2];

  double hxi;
  double hzt;

  switch (fce_nrm(f, v)) {
    case 0:
      hxi = v[f->vtx[1]]->y - v[f->vtx[0]]->y;
      hzt = v[f->vtx[2]]->z - v[f->vtx[0]]->z;

      break;
    case 1:
      hxi = v[f->vtx[1]]->x - v[f->vtx[0]]->x;
      hzt = v[f->vtx[2]]->z - v[f->vtx[0]]->z;

      break;
    case 2:
      hxi = v[f->vtx[1]]->x - v[f->vtx[0]]->x;
      hzt = v[f->vtx[2]]->y - v[f->vtx[0]]->y;

      break;
  }

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      mxi[i][j] = hxi * M[i][j];
      mzt[i][j] = hzt * M[i][j];
    }

  double* bloc;
  double dec[4];

  switch (f->cnd.type) {
    case DIR:
      break;
    case NEU:
      bloc = f->loc.neu.b->vp;

      for (int i = 0; i < 4; ++i)
        dec[i] = f->cnd.pps.neu.tta(v[f->vtx[i]]);

      for (int i = 0; i < 4; ++i) {
        bloc[i] = 0;

        int mui = MU[i];
        int nui = NU[i];

        for (int j = 0; j < 4; ++j) {
          int muj = MU[j];
          int nuj = NU[j];

          bloc[i] += dec[j] * mxi[mui][muj] * mzt[nui][nuj];
        }
      }

      break;
    case ROB:
      bloc = f->loc.rob.b->vp;

      double** mloc = f->loc.rob.m->v;
      double bet = f->cnd.pps.rob.bet;

      for (int i = 0; i < 4; ++i)
        dec[i] = f->cnd.pps.rob.tmp(v[f->vtx[i]]);

      for (int i = 0; i < 4; ++i) {
        bloc[i] = 0;

        int mui = MU[i];
        int nui = NU[i];

        for (int j = 0; j < 4; ++j) {
          int muj = MU[j];
          int nuj = NU[j];

          mloc[i][j] = bet * mxi[mui][muj] * mzt[nui][nuj];
          bloc[i] += dec[j] * mxi[mui][muj] * mzt[nui][nuj];
        }

        bloc[i] *= bet;
      }
  }

  return 0;
}

int fce_mov(struct fce* f, struct mtx_csj* a, struct vec* b) {
  double* bloc;

  switch (f->cnd.type) {
    case DIR:
      return 0;
    case NEU:
      bloc = f->loc.neu.b->vp;

      for (int i = 0; i < 4; ++i)
        b->vp[f->vtx[i]] += bloc[i];

      break;
    case ROB:
      bloc = f->loc.rob.b->vp;

      double** mloc = f->loc.rob.m->v;

      for (int i = 0; i < 4; ++i) {
        int iglob = f->vtx[i];

        b->vp[iglob] += bloc[i];

        for (int j = 0; j < 4; ++j) {
          int jglob = f->vtx[j];

          if (iglob == jglob) {
            a->dr[iglob] += mloc[i][j];
            continue;
          }

          if (iglob > jglob) {
            int lr0 = a->ia[iglob];
            int lr1 = a->ia[iglob + 1];

            for (int lr = lr0; lr < lr1; ++lr) {
              int lj = a->ja[lr];

              if (lj == jglob) {
                a->lr[lr] += mloc[i][j];
                break;
              }
            }
          } else {
            int ur0 = a->ia[jglob];
            int ur1 = a->ia[jglob + 1];

            for (int ur = ur0; ur < ur1; ++ur) {
              int ui = a->ja[ur];

              if (ui == iglob) {
                a->ur[ur] += mloc[i][j];
                break;
              }
            }
          }
        }
      }
  }

  return 0;
}
