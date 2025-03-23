#ifndef NUMX_COM_LOG_H
#define NUMX_COM_LOG_H

#include <numx/com/dev/log.h>

log_def(ilog, irec, int);
log_def(ulog, urec, unsigned);
log_def(flog, frec, float);
log_def(dlog, drec, double);
log_def(plog, prec, void *);

#define log_new(X)                                                                                                     \
    _Generic((X),                                                                                                      \
        struct ilog *: ilog_new,                                                                                       \
        struct ulog *: ulog_new,                                                                                       \
        struct flog *: flog_new,                                                                                       \
        struct dlog *: dlog_new,                                                                                       \
        struct plog *: plog_new)(X)

#define log_cls(X)                                                                                                     \
    _Generic((X),                                                                                                      \
        struct ilog *: ilog_cls,                                                                                       \
        struct ulog *: ulog_cls,                                                                                       \
        struct flog *: flog_cls,                                                                                       \
        struct dlog *: dlog_cls,                                                                                       \
        struct plog *: plog_cls)(X)

#define log_add(X, e)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct ilog *: ilog_add,                                                                                       \
        struct ulog *: ulog_add,                                                                                       \
        struct flog *: flog_add,                                                                                       \
        struct dlog *: dlog_add,                                                                                       \
        struct plog *: plog_add)(X, e)

#define log_ins(X, e, d)                                                                                               \
    _Generic((X),                                                                                                      \
        struct ilog *: ilog_ins,                                                                                       \
        struct ulog *: ulog_ins,                                                                                       \
        struct flog *: flog_ins,                                                                                       \
        struct dlog *: dlog_ins,                                                                                       \
        struct plog *: plog_ins)(X, e, d)

#define log_pop(X, e)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct ilog *: ilog_pop,                                                                                       \
        struct ulog *: ulog_pop,                                                                                       \
        struct flog *: flog_pop,                                                                                       \
        struct dlog *: dlog_pop,                                                                                       \
        struct plog *: plog_pop)(X, e)

#define log_beg(X, e)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct ilog *: ilog_beg,                                                                                       \
        struct ulog *: ulog_beg,                                                                                       \
        struct flog *: flog_beg,                                                                                       \
        struct dlog *: dlog_beg,                                                                                       \
        struct plog *: plog_beg)(X, e)

#define log_end(X, e)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct ilog *: ilog_end,                                                                                       \
        struct ulog *: ulog_end,                                                                                       \
        struct flog *: flog_end,                                                                                       \
        struct dlog *: dlog_end,                                                                                       \
        struct plog *: plog_end)(X, e)

#define log_rst(X)                                                                                                     \
    _Generic((X),                                                                                                      \
        struct ilog *: ilog_rst,                                                                                       \
        struct ulog *: ulog_rst,                                                                                       \
        struct flog *: flog_rst,                                                                                       \
        struct dlog *: dlog_rst,                                                                                       \
        struct plog *: plog_rst)(X)

#define log_adv(X, e)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct ilog *: ilog_adv,                                                                                       \
        struct ulog *: ulog_adv,                                                                                       \
        struct flog *: flog_adv,                                                                                       \
        struct dlog *: dlog_adv,                                                                                       \
        struct plog *: plog_adv)(X, e)

#define log_dec(X, e)                                                                                                  \
    _Generic((X),                                                                                                      \
        struct ilog *: ilog_dec,                                                                                       \
        struct ulog *: ulog_dec,                                                                                       \
        struct flog *: flog_dec,                                                                                       \
        struct dlog *: dlog_dec,                                                                                       \
        struct plog *: plog_dec)(X, e)

#endif // NUMX_COM_LOG_H
