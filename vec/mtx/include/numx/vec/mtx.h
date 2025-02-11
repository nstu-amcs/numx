#ifndef NUMX_VEC_MTX_H
#define NUMX_VEC_MTX_H

#include <numx/vec/vec.h>

struct imtx {
  struct imtx_pps {
    int m;
    int n;
  } pps;

  double** dat;
};

int imtx_new(struct imtx* m, struct imtx_pps pps);
int imtx_cls(struct imtx* m);

int imtx_vmlt(struct imtx* m, struct vec* x, struct vec* f);
int imtx_mmlt(struct imtx* a, struct imtx* b, struct imtx* r);

struct smtx {
  struct smtx_pps {
    int n;
    int z;
  } pps;

  double* dr;
  double* lr;
  double* ur;

  int* ia;
  int* ja;
};

int smtx_new(struct smtx* m, struct smtx_pps pps);
int smtx_cls(struct smtx* m);

int smtx_ilu(struct smtx* m, struct smtx* r);
int smtx_dgl(struct smtx* m, struct smtx* r);

int smtx_vmlt(struct smtx* m, struct vec* x, struct vec* f);
int smtx_mmlt(struct smtx* a, struct smtx* b, struct smtx* r);

// See https://www.ibm.com/docs/en/essl/6.2?topic=representation-compressed-diagonal-storage-mode
struct dmtx {
  struct dmtx_pps {
    int n;
    int d;
  } pps;

  double** ad;

  int* la;
};

int dmtx_new(struct dmtx* m, struct dmtx_pps pps);
int dmtx_cls(struct dmtx* m);

int dmtx_vmlt(struct dmtx* m, struct vec* x, struct vec* f);

struct jmtx {
  struct jmtx_pps {
    int m;
    int n;
  } pps;

  double (***dat)(struct vec*);
};

int jmtx_new(struct jmtx* m, struct jmtx_pps pps);
int jmtx_cls(struct jmtx* m);

#define mtx_new(X, p) _Generic((X), \
    struct imtx*: imtx_new,         \
    struct smtx*: smtx_new,         \
    struct dmtx*: dmtx_new,         \
    struct jmtx*: jmtx_new          \
    )(X, p)

#define mtx_cls(X) _Generic((X),  \
    struct imtx*: imtx_cls,       \
    struct smtx*: smtx_cls,       \
    struct dmtx*: dmtx_cls,       \
    struct jmtx*: jmtx_cls        \
    )(X)

#define mtx_ilu(X, r) _Generic((X), \
    struct smtx*: smtx_ilu          \
    )(X, r)

#define mtx_dgl(X, r) _Generic((X), \
    struct smtx*: smtx_dgl          \
    )(X, r)

#define mtx_vmlt(X, v, r) _Generic((X), \
    struct imtx*: imtx_vmlt,            \
    struct smtx*: smtx_vmlt,            \
    struct dmtx*: dmtx_vmlt             \
    )(X, v, r)

#define mtx_mmlt(X, b, r) _Generic((X), \
    struct imtx*: imtx_mmlt,            \
    struct smtx*: smtx_mmlt             \
    )(X, b, r)

#endif  // NUMX_VEC_MTX_H
