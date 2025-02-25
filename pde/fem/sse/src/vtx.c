#include <fem/sse/vtx.h>

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

int vtx_ini(struct vtx** h) {
  if (!h) {
    errno = EINVAL;
    return -1;
  }

  struct vtx* v = malloc(sizeof(struct vtx));

  if (!v) {
    errno = ENOMEM;
    return -1;
  }

  v->x = 0;
  v->y = 0;
  v->z = 0;

  *h = v;

  return 0;
}

int vtx_cls(struct vtx** h) {
  if (!h || !(*h)) {
    errno = EINVAL;
    return -1;
  }

  free(*h);
  *h = 0;

  return 0;
}

int vtx_sget(struct vtx* v, const char* buf) {
  if (buf[0] != 'v') {
    errno = EINVAL;
    return -1;
  }

  sscanf(buf + 1, "%lf %lf %lf", &v->x, &v->y, &v->z);

  return 0;
}
