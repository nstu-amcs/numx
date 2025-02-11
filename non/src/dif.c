#include <errno.h>
#include <numx/non/dif.h>

int pdif(double (*f)(struct vec*), int i, double hop, struct vec* p, double* r) {
  if (!f || !p || !r || i >= p->n) {
    errno = EINVAL;
    return -1;
  }

  p->dat[i] += hop;
  double f1 = f(p);

  p->dat[i] -= 2 * hop;
  double f2 = f(p);

  p->dat[i] += hop;
  *r = (f1 - f2) / (2 * hop);

  return 0;
}
