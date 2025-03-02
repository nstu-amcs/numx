#include <errno.h>
#include <numx/pde/cnd.h>
#include <stdio.h>

int cnd_get(struct cnd* cnd, const char* buf, struct pcut* dat) {
  if (!cnd || !buf) {
    errno = EINVAL;
    return -1;
  }

  char t = 0;
  int f = 0;
  int n = 0;

  if (sscanf(buf, "%c %d %n", &t, &f, &n) != 2)
    return -1;

  if (f >= dat->len) {
    errno = EINVAL;
    return -1;
  }

  switch (t) {
    case 'D':
      cnd->type = DIR;
      cnd->pps.dir.tmp = dat->dat[f];

      break;
    case 'N':
      cnd->type = NEU;
      cnd->pps.neu.tta = dat->dat[f];

      break;
    case 'R':
      cnd->type = ROB;
      cnd->pps.rob.tmp = dat->dat[f];

      if (sscanf(buf + n, "%lf", &cnd->pps.rob.bet) != 1)
        return -1;

      break;
    default:
      errno = EINVAL;
      return -1;
  }

  return 0;
}