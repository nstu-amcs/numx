#ifndef NUMX_COM_DEV_CUT_H
#define NUMX_COM_DEV_CUT_H

#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <numx/com/cap.h>

#define PUB
#define PRV static

#define cut_def(N, T)                                                                              \
    typedef struct N                                                                               \
    {                                                                                              \
        int cap;                                                                                   \
        int len;                                                                                   \
                                                                                                   \
        bool own;                                                                                  \
                                                                                                   \
        struct icap cmp;                                                                           \
        struct vcap anc;                                                                           \
                                                                                                   \
        T *dat;                                                                                    \
    } N;                                                                                           \
                                                                                                   \
    int N##_new(struct N *c);                                                                      \
    int N##_cls(struct N *c);                                                                      \
                                                                                                   \
    int N##_cov(struct N *c, T *s, int len);                                                       \
    int N##_mov(struct N *c, T *s, int len);                                                       \
                                                                                                   \
    int N##_shr(struct N *c);                                                                      \
    int N##_exp(struct N *c, int cap);                                                             \
    int N##_dev(struct N *c, int len);                                                             \
                                                                                                   \
    int N##_add(struct N *c, T e);                                                                 \
    int N##_srt(struct N *c);                                                                      \
                                                                                                   \
    extern int _unused_

#define cut_gen(N, T, M)                                                                           \
    M int N##_new(struct N *c)                                                                     \
    {                                                                                              \
        c->cap = 0;                                                                                \
        c->len = 0;                                                                                \
        c->own = 0;                                                                                \
        c->dat = 0;                                                                                \
        c->cmp.ctx = 0;                                                                            \
        c->anc.ctx = 0;                                                                            \
        c->cmp.run = 0;                                                                            \
        c->anc.run = 0;                                                                            \
                                                                                                   \
        return 0;                                                                                  \
    }                                                                                              \
                                                                                                   \
    M int N##_cls(struct N *c)                                                                     \
    {                                                                                              \
        if (!c) {                                                                                  \
            errno = EINVAL;                                                                        \
            return -1;                                                                             \
        }                                                                                          \
                                                                                                   \
        if (c->dat && c->own) {                                                                    \
            free(c->dat);                                                                          \
            c->dat = 0;                                                                            \
            c->cap = 0;                                                                            \
            c->len = 0;                                                                            \
        }                                                                                          \
                                                                                                   \
        return 0;                                                                                  \
    }                                                                                              \
                                                                                                   \
    M int N##_cov(struct N *c, T *s, int len)                                                      \
    {                                                                                              \
        if (!c || !s || len == 0) {                                                                \
            errno = EINVAL;                                                                        \
            return -1;                                                                             \
        }                                                                                          \
                                                                                                   \
        if (c->dat && c->own)                                                                      \
            free(c->dat);                                                                          \
                                                                                                   \
        c->own = false;                                                                            \
        c->cap = len;                                                                              \
        c->len = len;                                                                              \
        c->dat = s;                                                                                \
                                                                                                   \
        return 0;                                                                                  \
    }                                                                                              \
                                                                                                   \
    M int N##_mov(struct N *c, T *s, int len)                                                      \
    {                                                                                              \
        if (!c || !s) {                                                                            \
            errno = EINVAL;                                                                        \
            return -1;                                                                             \
        }                                                                                          \
                                                                                                   \
        if (c->dat && c->own)                                                                      \
            free(c->dat);                                                                          \
                                                                                                   \
        c->own = true;                                                                             \
        c->cap = len;                                                                              \
        c->len = len;                                                                              \
        c->dat = s;                                                                                \
                                                                                                   \
        return 0;                                                                                  \
    }                                                                                              \
                                                                                                   \
    M int N##_shr(struct N *c)                                                                     \
    {                                                                                              \
        if (!c) {                                                                                  \
            errno = EINVAL;                                                                        \
            return -1;                                                                             \
        }                                                                                          \
                                                                                                   \
        if (!c->own || c->len == c->cap)                                                           \
            return 0;                                                                              \
                                                                                                   \
        if (c->len == 0) {                                                                         \
            free(c->dat);                                                                          \
                                                                                                   \
            c->cap = 0;                                                                            \
            c->own = false;                                                                        \
            c->dat = 0;                                                                            \
                                                                                                   \
            return 0;                                                                              \
        }                                                                                          \
                                                                                                   \
        T *n = realloc(c->dat, sizeof(T) * c->len);                                                \
                                                                                                   \
        if (!n)                                                                                    \
            return -1;                                                                             \
                                                                                                   \
        c->cap = c->len;                                                                           \
        c->dat = n;                                                                                \
                                                                                                   \
        return 0;                                                                                  \
    }                                                                                              \
                                                                                                   \
    M int N##_exp(struct N *c, int cap)                                                            \
    {                                                                                              \
        if (cap <= c->cap)                                                                         \
            return 0;                                                                              \
                                                                                                   \
        if (!c->dat || !c->own) {                                                                  \
            T *n = malloc(sizeof(T) * cap);                                                        \
                                                                                                   \
            if (!n) {                                                                              \
                errno = ENOMEM;                                                                    \
                return -1;                                                                         \
            }                                                                                      \
                                                                                                   \
            c->cap = cap;                                                                          \
                                                                                                   \
            if (c->dat)                                                                            \
                memcpy(n, c->dat, sizeof(T) * c->len);                                             \
                                                                                                   \
            c->dat = n;                                                                            \
            c->own = true;                                                                         \
                                                                                                   \
            return 0;                                                                              \
        }                                                                                          \
                                                                                                   \
        T *n = realloc(c->dat, sizeof(T) * cap);                                                   \
                                                                                                   \
        if (!n) {                                                                                  \
            errno = ENOMEM;                                                                        \
            return -1;                                                                             \
        }                                                                                          \
                                                                                                   \
        c->cap = cap;                                                                              \
        c->dat = n;                                                                                \
                                                                                                   \
        return 0;                                                                                  \
    }                                                                                              \
                                                                                                   \
    M int N##_dev(struct N *c, int num)                                                            \
    {                                                                                              \
        if (!c || num < 1) {                                                                       \
            errno = EINVAL;                                                                        \
            return -1;                                                                             \
        }                                                                                          \
                                                                                                   \
        if (c->len + num > c->cap)                                                                 \
            if (N##_exp(c, c->len + num))                                                          \
                return -1;                                                                         \
                                                                                                   \
        memset(c->dat + c->len, 0, sizeof(T) * num);                                               \
        c->len += num;                                                                             \
                                                                                                   \
        return 0;                                                                                  \
    }                                                                                              \
                                                                                                   \
    M int N##_add(struct N *c, T e)                                                                \
    {                                                                                              \
        if (!c) {                                                                                  \
            errno = EINVAL;                                                                        \
            return -1;                                                                             \
        }                                                                                          \
                                                                                                   \
        if (c->len == c->cap)                                                                      \
            if (N##_exp(c, (c->cap + 1) * 2))                                                      \
                return -1;                                                                         \
                                                                                                   \
        c->dat[c->len++] = e;                                                                      \
                                                                                                   \
        return 0;                                                                                  \
    }                                                                                              \
                                                                                                   \
    M int N##_srt(struct N *c)                                                                     \
    {                                                                                              \
        if (!c) {                                                                                  \
            errno = EINVAL;                                                                        \
            return -1;                                                                             \
        }                                                                                          \
                                                                                                   \
        errno = ENOTSUP;                                                                           \
                                                                                                   \
        return -1;                                                                                 \
    }                                                                                              \
                                                                                                   \
    extern int _unused_

#endif // NUMX_COM_DEV_CUT_H
