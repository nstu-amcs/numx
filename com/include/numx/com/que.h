#ifndef NUMX_COM_QUE_H
#define NUMX_COM_QUE_H

#include <numx/com/cut.h>
#include <numx/com/dev/que.h>

que_def(ique, icut, int);
que_def(uque, ucut, unsigned);
que_def(fque, fcut, float);
que_def(dque, dcut, double);
que_def(pque, pcut, void *);

#define que_fixu(X, i)                                                                                                 \
    _Generic((X),                                                                                                      \
        struct icut *: ique_fixu,                                                                                      \
        struct ucut *: uque_fixu,                                                                                      \
        struct fcut *: fque_fixu,                                                                                      \
        struct dcut *: dque_fixu,                                                                                      \
        struct pcut *: pque_fixu)(X, i)

#define que_fixd(X, i)                                                                                                 \
    _Generic((X),                                                                                                      \
        struct icut *: ique_fixd,                                                                                      \
        struct ucut *: uque_fixd,                                                                                      \
        struct fcut *: fque_fixd,                                                                                      \
        struct dcut *: dque_fixd,                                                                                      \
        struct pcut *: pque_fixd)(X, i)

#define que_fix(X)                                                                                                     \
    _Generic((X),                                                                                                      \
        struct icut *: ique_fix,                                                                                       \
        struct ucut *: uque_fix,                                                                                       \
        struct fcut *: fque_fix,                                                                                       \
        struct dcut *: dque_fix,                                                                                       \
        struct pcut *: pque_fix)(X)

#define que_add(X, e)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct icut *: ique_add,                                                                                       \
        struct ucut *: uque_add,                                                                                       \
        struct fcut *: fque_add,                                                                                       \
        struct dcut *: dque_add,                                                                                       \
        struct pcut *: pque_add)(X, e)

#define que_ext(X, e)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct icut *: ique_ext,                                                                                       \
        struct ucut *: uque_ext,                                                                                       \
        struct fcut *: fque_ext,                                                                                       \
        struct dcut *: dque_ext,                                                                                       \
        struct pcut *: pque_ext)(X, e)

#define que_srt(X)                                                                                                     \
    _Generic((X),                                                                                                      \
        struct icut *: ique_srt,                                                                                       \
        struct ucut *: uque_srt,                                                                                       \
        struct fcut *: fque_srt,                                                                                       \
        struct dcut *: dque_srt,                                                                                       \
        struct pcut *: pque_srt)(X)

#endif // NUMX_COM_QUE_H
