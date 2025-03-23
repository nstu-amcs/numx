#ifndef NUMX_COM_CUT_H
#define NUMX_COM_CUT_H

#include <numx/com/dev/cut.h>

cut_def(icut, int);
cut_def(ucut, unsigned);
cut_def(fcut, float);
cut_def(dcut, double);
cut_def(pcut, void *);

#define cut_new(X)                                                                                                     \
    _Generic((X),                                                                                                      \
        struct icut *: icut_new,                                                                                       \
        struct ucut *: ucut_new,                                                                                       \
        struct fcut *: fcut_new,                                                                                       \
        struct dcut *: dcut_new,                                                                                       \
        struct pcut *: pcut_new)(X)

#define cut_cls(X)                                                                                                     \
    _Generic((X),                                                                                                      \
        struct icut *: icut_cls,                                                                                       \
        struct ucut *: ucut_cls,                                                                                       \
        struct fcut *: fcut_cls,                                                                                       \
        struct dcut *: dcut_cls,                                                                                       \
        struct pcut *: pcut_cls)(X)

#define cut_cov(X, s, l)                                                                                               \
    _Generic((X),                                                                                                      \
        struct icut *: icut_cov,                                                                                       \
        struct ucut *: ucut_cov,                                                                                       \
        struct fcut *: fcut_cov,                                                                                       \
        struct dcut *: dcut_cov,                                                                                       \
        struct pcut *: pcut_cov)(X, s, l)

#define cut_mov(X, s, l)                                                                                               \
    _Generic((X),                                                                                                      \
        struct icut *: icut_mov,                                                                                       \
        struct ucut *: ucut_mov,                                                                                       \
        struct fcut *: fcut_mov,                                                                                       \
        struct dcut *: dcut_mov,                                                                                       \
        struct pcut *: pcut_mov)(X, s, l)

#define cut_shr(X)                                                                                                     \
    _Generic((X),                                                                                                      \
        struct icut *: icut_shr,                                                                                       \
        struct ucut *: ucut_shr,                                                                                       \
        struct fcut *: fcut_shr,                                                                                       \
        struct dcut *: dcut_shr,                                                                                       \
        struct pcut *: pcut_shr)(X)

#define cut_exp(X, c)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct icut *: icut_exp,                                                                                       \
        struct ucut *: ucut_exp,                                                                                       \
        struct fcut *: fcut_exp,                                                                                       \
        struct dcut *: dcut_exp,                                                                                       \
        struct pcut *: pcut_exp)(X, c)

#define cut_dev(X, l)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct icut *: icut_dev,                                                                                       \
        struct ucut *: ucut_dev,                                                                                       \
        struct fcut *: fcut_dev,                                                                                       \
        struct dcut *: dcut_dev,                                                                                       \
        struct pcut *: pcut_dev)(X, l)

#define cut_add(X, e)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct icut *: icut_add,                                                                                       \
        struct ucut *: ucut_add,                                                                                       \
        struct fcut *: fcut_add,                                                                                       \
        struct dcut *: dcut_add,                                                                                       \
        struct pcut *: pcut_add)(X, e)

#define cut_srt(X)                                                                                                     \
    _Generic((X),                                                                                                      \
        struct icut *: icut_srt,                                                                                       \
        struct ucut *: ucut_srt,                                                                                       \
        struct fcut *: fcut_srt,                                                                                       \
        struct dcut *: dcut_srt,                                                                                       \
        struct pcut *: pcut_srt)(X)

#endif // NUMX_COM_CUT_H
