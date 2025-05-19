#ifndef NUMX_VEC_MTX_H
#define NUMX_VEC_MTX_H

#include <numx/non/fun.h>
#include <numx/vec/vec.h>

/** Matrix in dense storage mode. */
typedef struct imtx
{
    struct imtx_pps
    {
        int m;
        int n;
    } pps;

    double **dat;
} imtx;

int imtx_new(struct imtx *m, struct imtx_pps pps);
int imtx_cls(struct imtx *m);

int imtx_inc(struct imtx *m, int i, int j, double v);
int imtx_cmb(struct imtx *a, struct imtx *b, struct imtx *r, double k);

int imtx_vdup(struct imtx *s, struct imtx *d);
int imtx_sdup(struct imtx *s, struct imtx *d);

int imtx_rst(struct imtx *m);

int imtx_ilu(struct imtx *m, struct imtx *r);
int imtx_dgl(struct imtx *m, struct imtx *r);

int imtx_vmul(struct imtx *m, struct vec *x, struct vec *f);
int imtx_mmul(struct imtx *a, struct imtx *b, struct imtx *r);

/**
 * Sparse matrix in row-column storage mode.
 *
 * See https://www.ibm.com/docs/en/essl/6.2?topic=representation-storage-by-rows
 */
typedef struct smtx
{
    struct smtx_pps
    {
        int n;
        int z;
    } pps;

    double *dr;
    double *lr;
    double *ur;

    int *ia;
    int *ja;
} smtx;

int smtx_new(struct smtx *m, struct smtx_pps pps);
int smtx_cls(struct smtx *m);

int smtx_inc(struct smtx *m, int i, int j, double v);
int smtx_cmb(struct smtx *a, struct smtx *b, struct smtx *r, double k);

int smtx_vdup(struct smtx *s, struct smtx *d);
int smtx_sdup(struct smtx *s, struct smtx *d);

int smtx_rst(struct smtx *m);

int smtx_ilu(struct smtx *m, struct smtx *r);
int smtx_dgl(struct smtx *m, struct smtx *r);

int smtx_vmul(struct smtx *m, struct vec *x, struct vec *f);
int smtx_mmul(struct smtx *a, struct smtx *b, struct smtx *r);

/**
 * Sparse matrix in compressed-diagonal storage mode.
 *
 * See
 * https://www.ibm.com/docs/en/essl/6.2?topic=representation-compressed-diagonal-storage-mode
 * */
typedef struct dmtx
{
    struct dmtx_pps
    {
        int n;
        int d;
    } pps;

    double **ad;

    int *la;
} dmtx;

int dmtx_new(struct dmtx *m, struct dmtx_pps pps);
int dmtx_cls(struct dmtx *m);

int dmtx_inc(struct dmtx *m, int i, int j, double v);
int dmtx_cmb(struct dmtx *a, struct dmtx *b, struct dmtx *r, double k);

int dmtx_vdup(struct dmtx *s, struct dmtx *d);
int dmtx_sdup(struct dmtx *s, struct dmtx *d);

int dmtx_rst(struct dmtx *m);

int dmtx_ilu(struct dmtx *m, struct dmtx *r);
int dmtx_dgl(struct dmtx *m, struct dmtx *r);

int dmtx_vmul(struct dmtx *m, struct vec *x, struct vec *f);
int dmtx_mmul(struct dmtx *a, struct dmtx *b, struct dmtx *r);

typedef struct jmtx
{
    struct jmtx_pps
    {
        int m;
        int n;
    } pps;

    mfun **dat;
} jmtx;

int jmtx_new(struct jmtx *m, struct jmtx_pps pps);
int jmtx_cls(struct jmtx *m);

#define mtx_new(X, p)                                                          \
    _Generic((X),                                                              \
        struct imtx *: imtx_new,                                               \
        struct smtx *: smtx_new,                                               \
        struct dmtx *: dmtx_new,                                               \
        struct jmtx *: jmtx_new)(X, p)

#define mtx_cls(X)                                                             \
    _Generic((X),                                                              \
        struct imtx *: imtx_cls,                                               \
        struct smtx *: smtx_cls,                                               \
        struct dmtx *: dmtx_cls,                                               \
        struct jmtx *: jmtx_cls)(X)

#define mtx_inc(X, i, j, v)                                                    \
    _Generic((X),                                                              \
        struct imtx *: imtx_inc,                                               \
        struct smtx *: smtx_inc,                                               \
        struct dmtx *: dmtx_inc)(X, i, j, v)

#define mtx_cmb(X, b, r, k)                                                    \
    _Generic((X),                                                              \
        struct imtx *: imtx_cmb,                                               \
        struct smtx *: smtx_cmb,                                               \
        struct dmtx *: dmtx_cmb)(X, b, r, k)

#define mtx_ilu(X, r)                                                          \
    _Generic((X),                                                              \
        struct imtx *: imtx_ilu,                                               \
        struct smtx *: smtx_ilu,                                               \
        struct dmtx *: dmtx_ilu)(X, r)

#define mtx_dgl(X, r)                                                          \
    _Generic((X),                                                              \
        struct imtx *: imtx_dgl,                                               \
        struct smtx *: smtx_dgl,                                               \
        struct dmtx *: dmtx_dgl)(X, r)

#define mtx_vdup(X, d)                                                         \
    _Generic((X),                                                              \
        struct imtx *: imtx_vdup,                                              \
        struct smtx *: smtx_vdup,                                              \
        struct dmtx *: dmtx_vdup)(X, d)

#define mtx_sdup(X, d)                                                         \
    _Generic((X),                                                              \
        struct imtx *: imtx_sdup,                                              \
        struct smtx *: smtx_sdup,                                              \
        struct dmtx *: dmtx_sdup)(X, d)

#define mtx_rst(X)                                                             \
    _Generic((X),                                                              \
        struct imtx *: imtx_rst,                                               \
        struct smtx *: smtx_rst,                                               \
        struct dmtx *: dmtx_rst)(X)

#define mtx_vmul(X, v, r)                                                      \
    _Generic((X),                                                              \
        struct imtx *: imtx_vmul,                                              \
        struct smtx *: smtx_vmul,                                              \
        struct dmtx *: dmtx_vmul)(X, v, r)

#define mtx_mmul(X, m, r)                                                      \
    _Generic((X),                                                              \
        struct imtx *: imtx_mmul,                                              \
        struct smtx *: smtx_mmul,                                              \
        struct dmtx *: dmtx_mmul)(X, m, r)

#endif // NUMX_VEC_MTX_H
