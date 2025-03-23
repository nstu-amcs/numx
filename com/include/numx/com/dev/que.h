#ifndef NUMX_COM_DEV_QUE_H
#define NUMX_COM_DEV_QUE_H

#define PQUE_L(i) (i * 2 + 1)
#define PQUE_R(i) (i * 2 + 2)
#define PQUE_P(i) ((i - 1) / 2)

#define que_def(N, C, T)                                                                                               \
    int N##_fixu(struct C *q, int i);                                                                                  \
    int N##_fixd(struct C *q, int i);                                                                                  \
                                                                                                                       \
    int N##_fix(struct C *q);                                                                                          \
    int N##_add(struct C *q, T e);                                                                                     \
    int N##_ext(struct C *q, T *e);                                                                                    \
    int N##_srt(struct C *q);                                                                                          \
                                                                                                                       \
    extern int _unused_

#define que_gen(N, C, T, M)                                                                                            \
    static inline void N##_swap(struct C *c, int a, int b)                                                             \
    {                                                                                                                  \
        T t = c->dat[a];                                                                                               \
                                                                                                                       \
        c->dat[a] = c->dat[b];                                                                                         \
        c->dat[b] = t;                                                                                                 \
    }                                                                                                                  \
                                                                                                                       \
    M int N##_fixu(struct C *q, int i)                                                                                 \
    {                                                                                                                  \
        if (!q) {                                                                                                      \
            errno = EINVAL;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        int p = PQUE_P(i);                                                                                             \
                                                                                                                       \
        while (i > 0 && q->cmp.run(q->cmp.ctx, 2, q->dat[i], q->dat[p]) == -1) {                                       \
            N##_swap(q, p, i);                                                                                         \
                                                                                                                       \
            i = p;                                                                                                     \
            p = PQUE_P(i);                                                                                             \
        }                                                                                                              \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int N##_fixd(struct C *q, int i)                                                                                 \
    {                                                                                                                  \
        if (!q) {                                                                                                      \
            errno = EINVAL;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        if (i >= q->len) {                                                                                             \
            errno = ERANGE;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        int l = PQUE_L(i);                                                                                             \
        int r = PQUE_R(i);                                                                                             \
        int m;                                                                                                         \
                                                                                                                       \
        if (l < q->len && q->cmp.run(q->cmp.ctx, 2, q->dat[l], q->dat[i]) == -1)                                       \
            m = l;                                                                                                     \
        else                                                                                                           \
            m = i;                                                                                                     \
                                                                                                                       \
        if (r < q->len && q->cmp.run(q->cmp.ctx, 2, q->dat[r], q->dat[m]) == -1)                                       \
            m = r;                                                                                                     \
                                                                                                                       \
        if (m != i) {                                                                                                  \
            N##_swap(q, m, i);                                                                                         \
                                                                                                                       \
            return N##_fixd(q, m);                                                                                     \
        }                                                                                                              \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int N##_fix(struct C *q)                                                                                         \
    {                                                                                                                  \
        if (!q) {                                                                                                      \
            errno = EINVAL;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        for (int i = q->len / 2 - 1;; --i) {                                                                           \
            if (N##_fixd(q, i))                                                                                        \
                return -1;                                                                                             \
                                                                                                                       \
            if (i == 0)                                                                                                \
                break;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int N##_add(struct C *q, T e)                                                                                    \
    {                                                                                                                  \
        if (!q) {                                                                                                      \
            errno = EINVAL;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        if (C##_add(q, e))                                                                                             \
            return -1;                                                                                                 \
                                                                                                                       \
        return N##_fixu(q, q->len - 1);                                                                                \
    }                                                                                                                  \
                                                                                                                       \
    M int N##_ext(struct C *q, T *e)                                                                                   \
    {                                                                                                                  \
        if (!q) {                                                                                                      \
            errno = EINVAL;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        if (q->len == 0) {                                                                                             \
            errno = ENOENT;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        *e = q->dat[0];                                                                                                \
        N##_swap(q, 0, q->len - 1);                                                                                    \
        q->len -= 1;                                                                                                   \
                                                                                                                       \
        if (q->len != 0)                                                                                               \
            return N##_fixd(q, 0);                                                                                     \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    M int N##_srt(struct C *q)                                                                                         \
    {                                                                                                                  \
        if (!q) {                                                                                                      \
            errno = EINVAL;                                                                                            \
            return -1;                                                                                                 \
        }                                                                                                              \
                                                                                                                       \
        if (q->len == 0)                                                                                               \
            return 0;                                                                                                  \
                                                                                                                       \
        int len = q->len;                                                                                              \
                                                                                                                       \
        for (int i = q->len - 1; i > 0; --i) {                                                                         \
            N##_swap(q, 0, i);                                                                                         \
            q->len -= 1;                                                                                               \
            N##_fixd(q, 0);                                                                                            \
        }                                                                                                              \
                                                                                                                       \
        q->len = len;                                                                                                  \
                                                                                                                       \
        return 0;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    extern int _unused_

#endif // NUMX_COM_DEV_QUE_H
