#ifndef NUMX_COM_DEV_LOG_H
#define NUMX_COM_DEV_LOG_H

#include <assert.h>
#include <errno.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>

#include <numx/com/cap.h>

#define PUB
#define PRV static

typedef enum log_ref
{
    LOG_L,
    LOG_I,
    LOG_R,
} log_ref;

#define log_def(L, R, T)                                                                                               \
    typedef struct R                                                                                                   \
    {                                                                                                                  \
        T e;                                                                                                           \
                                                                                                                       \
        struct R *next;                                                                                                \
        struct R *prev;                                                                                                \
    } R;                                                                                                               \
                                                                                                                       \
    typedef struct L                                                                                                   \
    {                                                                                                                  \
        int len;                                                                                                       \
                                                                                                                       \
        bool srt;                                                                                                      \
        bool dup;                                                                                                      \
                                                                                                                       \
        struct icap cmp;                                                                                               \
                                                                                                                       \
        struct R *beg;                                                                                                 \
        struct R *end;                                                                                                 \
        struct R *itr;                                                                                                 \
    } L;                                                                                                               \
                                                                                                                       \
    int L##_new(struct L *l);                                                                                          \
    int L##_cls(struct L *l);                                                                                          \
                                                                                                                       \
    int L##_add(struct L *l, T e);                                                                                     \
    int L##_ins(struct L *l, T e, enum log_ref r);                                                                     \
    int L##_pop(struct L *l, T *e);                                                                                    \
                                                                                                                       \
    int L##_beg(struct L *l, T *e);                                                                                    \
    int L##_end(struct L *l, T *e);                                                                                    \
                                                                                                                       \
    int L##_rst(struct L *l);                                                                                          \
    int L##_adv(struct L *l, T *e);                                                                                    \
    int L##_dec(struct L *l, T *e);                                                                                    \
                                                                                                                       \
    extern int _unused_

#define log_gen(L, R, T, M)                                                                                            \
    M int L##_new(struct L *l)                                                                                         \
    {                                                                                                                  \
        assert(l);                                                                                                     \
                                                                                                                       \
        l->len = 0;                                                                                                    \
        l->srt = false;                                                                                                \
        l->dup = true;                                                                                                 \
        l->cmp.run = 0;                                                                                                \
        l->cmp.ctx = 0;                                                                                                \
        l->beg = 0;                                                                                                    \
        l->end = 0;                                                                                                    \
        l->itr = 0;                                                                                                    \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int L##_cls(struct L *l)                                                                                         \
    {                                                                                                                  \
        assert(l);                                                                                                     \
                                                                                                                       \
        struct R *r = l->beg;                                                                                          \
        struct R *t;                                                                                                   \
                                                                                                                       \
        while (r) {                                                                                                    \
            t = r;                                                                                                     \
            r = r->next;                                                                                               \
                                                                                                                       \
            free(t);                                                                                                   \
        }                                                                                                              \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int L##_add(struct L *l, T e)                                                                                    \
    {                                                                                                                  \
        assert(l);                                                                                                     \
                                                                                                                       \
        l->itr = 0;                                                                                                    \
                                                                                                                       \
        if (l->srt || !l->dup)                                                                                         \
            for (T i; !L##_adv(l, &i);) {                                                                              \
                int c = l->cmp.run(l->cmp.ctx, 2, i, e);                                                               \
                                                                                                                       \
                if (c == 0 && !l->dup) {                                                                               \
                    errno = EALREADY;                                                                                  \
                    return -1;                                                                                         \
                }                                                                                                      \
                                                                                                                       \
                if (c != 1 && l->srt) {                                                                                \
                    if (L##_ins(l, e, LOG_L))                                                                          \
                        return -1;                                                                                     \
                                                                                                                       \
                    l->itr = l->itr->prev;                                                                             \
                                                                                                                       \
                    return 0;                                                                                          \
                }                                                                                                      \
            }                                                                                                          \
                                                                                                                       \
        if (L##_ins(l, e, LOG_R))                                                                                      \
            return -1;                                                                                                 \
                                                                                                                       \
        if (!l->itr)                                                                                                   \
            l->itr = l->end;                                                                                           \
        else                                                                                                           \
            l->itr = l->itr->next;                                                                                     \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int L##_ins(struct L *l, T e, enum log_ref d)                                                                    \
    {                                                                                                                  \
        assert(l);                                                                                                     \
                                                                                                                       \
        if (d == LOG_I) {                                                                                              \
            if (!l->itr) {                                                                                             \
                errno = ENOENT;                                                                                        \
                return -1;                                                                                             \
            }                                                                                                          \
                                                                                                                       \
            l->itr->e = e;                                                                                             \
                                                                                                                       \
            return 0;                                                                                                  \
        }                                                                                                              \
                                                                                                                       \
        struct R *r = malloc(sizeof(struct R));                                                                        \
                                                                                                                       \
        if (!r) {                                                                                                      \
            errno = ENOMEM;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        r->e = e;                                                                                                      \
                                                                                                                       \
        switch (d) {                                                                                                   \
            case LOG_I:                                                                                                \
            case LOG_L:                                                                                                \
                if (l->itr) {                                                                                          \
                    r->prev = l->itr->prev;                                                                            \
                    r->next = l->itr;                                                                                  \
                    l->itr->prev = r;                                                                                  \
                                                                                                                       \
                    if (!r->prev)                                                                                      \
                        l->beg = r;                                                                                    \
                    else                                                                                               \
                        r->prev->next = r;                                                                             \
                } else {                                                                                               \
                    r->prev = 0;                                                                                       \
                    r->next = l->beg;                                                                                  \
                    l->beg = r;                                                                                        \
                                                                                                                       \
                    if (!r->next)                                                                                      \
                        l->end = r;                                                                                    \
                    else                                                                                               \
                        r->next->prev = r;                                                                             \
                }                                                                                                      \
                                                                                                                       \
                break;                                                                                                 \
            case LOG_R:                                                                                                \
                if (l->itr) {                                                                                          \
                    r->next = l->itr->next;                                                                            \
                    r->prev = l->itr;                                                                                  \
                    l->itr->next = r;                                                                                  \
                                                                                                                       \
                    if (!r->next)                                                                                      \
                        l->end = r;                                                                                    \
                    else                                                                                               \
                        r->next->prev = r;                                                                             \
                } else {                                                                                               \
                    r->next = 0;                                                                                       \
                    r->prev = l->end;                                                                                  \
                    l->end = r;                                                                                        \
                                                                                                                       \
                    if (!r->prev)                                                                                      \
                        l->beg = r;                                                                                    \
                    else                                                                                               \
                        r->prev->next = r;                                                                             \
                }                                                                                                      \
                                                                                                                       \
                break;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        l->len += 1;                                                                                                   \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int L##_pop(struct L *l, T *e)                                                                                   \
    {                                                                                                                  \
        assert(l);                                                                                                     \
                                                                                                                       \
        struct R *r = l->itr;                                                                                          \
                                                                                                                       \
        if (!r) {                                                                                                      \
            errno = ENOENT;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        l->itr = r->prev;                                                                                              \
        l->len -= 1;                                                                                                   \
                                                                                                                       \
        if (e)                                                                                                         \
            *e = r->e;                                                                                                 \
                                                                                                                       \
        if (!r->prev)                                                                                                  \
            l->beg = r->next;                                                                                          \
        else                                                                                                           \
            r->prev->next = r->next;                                                                                   \
                                                                                                                       \
        if (!r->next)                                                                                                  \
            l->end = r->prev;                                                                                          \
        else                                                                                                           \
            r->next->prev = r->prev;                                                                                   \
                                                                                                                       \
        free(r);                                                                                                       \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int L##_beg(struct L *l, T *e)                                                                                   \
    {                                                                                                                  \
        assert(l);                                                                                                     \
                                                                                                                       \
        l->itr = l->beg;                                                                                               \
                                                                                                                       \
        if (!l->itr) {                                                                                                 \
            errno = ENOENT;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        if (e)                                                                                                         \
            *e = l->itr->e;                                                                                            \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int L##_end(struct L *l, T *e)                                                                                   \
    {                                                                                                                  \
        assert(l);                                                                                                     \
                                                                                                                       \
        l->itr = l->end;                                                                                               \
                                                                                                                       \
        if (!l->itr) {                                                                                                 \
            errno = ENOENT;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        if (e)                                                                                                         \
            *e = l->itr->e;                                                                                            \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int L##_rst(struct L *l)                                                                                         \
    {                                                                                                                  \
        assert(l);                                                                                                     \
                                                                                                                       \
        l->itr = 0;                                                                                                    \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int L##_adv(struct L *l, T *e)                                                                                   \
    {                                                                                                                  \
        assert(l);                                                                                                     \
                                                                                                                       \
        if (!l->itr) {                                                                                                 \
            l->itr = l->beg;                                                                                           \
                                                                                                                       \
            if (!l->itr) {                                                                                             \
                errno = ENOENT;                                                                                        \
                return -1;                                                                                             \
            }                                                                                                          \
        } else {                                                                                                       \
            if (!l->itr->next) {                                                                                       \
                errno = ENOENT;                                                                                        \
                return -1;                                                                                             \
            }                                                                                                          \
                                                                                                                       \
            l->itr = l->itr->next;                                                                                     \
        }                                                                                                              \
                                                                                                                       \
        if (e)                                                                                                         \
            *e = l->itr->e;                                                                                            \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int L##_dec(struct L *l, T *e)                                                                                   \
    {                                                                                                                  \
        assert(l);                                                                                                     \
                                                                                                                       \
        if (!l->itr) {                                                                                                 \
            l->itr = l->end;                                                                                           \
                                                                                                                       \
            if (!l->itr) {                                                                                             \
                errno = ENOENT;                                                                                        \
                return -1;                                                                                             \
            }                                                                                                          \
        } else {                                                                                                       \
            if (!l->itr->prev) {                                                                                       \
                errno = ENOENT;                                                                                        \
                return -1;                                                                                             \
            }                                                                                                          \
                                                                                                                       \
            l->itr = l->itr->prev;                                                                                     \
        }                                                                                                              \
                                                                                                                       \
        if (e)                                                                                                         \
            *e = l->itr->e;                                                                                            \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int L##_len(struct L *l)                                                                                         \
    {                                                                                                                  \
        assert(l);                                                                                                     \
        return l->len;                                                                                                 \
    }                                                                                                                  \
                                                                                                                       \
    extern int _unused_

#endif // NUMX_COM_DEV_LOG_H
