#include <math.h>
#include <stdio.h>

#include <numx/vec/dss.h>
#include <numx/vec/mtx.h>
#include <numx/vec/vec.h>

#define COND 0.1

double a1d[3] = {0, -500, 0};
double b1d[3] = {100, -500, 0};
double a2d[3] = {0, 0, 0};
double b2d[3] = {100, 0, 0};
double a3d[3] = {0, 500, 0};
double b3d[3] = {100, 500, 0};

double m1d[3] = {200, 0, 0};
double n1d[3] = {300, 0, 0};
double m2d[3] = {500, 0, 0};
double n2d[3] = {600, 0, 0};
double m3d[3] = {1000, 0, 0};
double n3d[3] = {1100, 0, 0};

struct vec a1 = {.n = 3, .dat = a1d};
struct vec b1 = {.n = 3, .dat = b1d};
struct vec a2 = {.n = 3, .dat = a2d};
struct vec b2 = {.n = 3, .dat = b2d};
struct vec a3 = {.n = 3, .dat = a3d};
struct vec b3 = {.n = 3, .dat = b3d};

struct vec m1 = {.n = 3, .dat = m1d};
struct vec n1 = {.n = 3, .dat = n1d};
struct vec m2 = {.n = 3, .dat = m2d};
struct vec n2 = {.n = 3, .dat = n2d};
struct vec m3 = {.n = 3, .dat = m3d};
struct vec n3 = {.n = 3, .dat = n3d};

struct lin
{
    struct vec *g1;
    struct vec *g2;
};

// Источник
struct sup
{
    struct lin l;
};

// Приемник
struct rec
{
    struct lin l;
    struct vec k;
};

struct sup ab1 = {
    .l = {.g1 = &a1, .g2 = &b1}
};
struct sup ab2 = {
    .l = {.g1 = &a2, .g2 = &b2}
};
struct sup ab3 = {
    .l = {.g1 = &a3, .g2 = &b3}
};

struct rec mn1 = {
    .l = {.g1 = &m1, .g2 = &n1}
};
struct rec mn2 = {
    .l = {.g1 = &m2, .g2 = &n2},
};
struct rec mn3 = {
    .l = {.g1 = &m3, .g2 = &n3}
};

struct sup *sup[] = {&ab1, &ab2, &ab3};
struct rec *rec[] = {&mn1, &mn2, &mn3};

// Расчет коэффициента установки
double k(struct lin *s, struct lin *r)
{
    double rbm = 0;
    double ram = 0;
    double rbn = 0;
    double ran = 0;

    // vec_dst - расстояние между векторами

    vec_dst(s->g2, r->g1, &rbm);
    vec_dst(s->g1, r->g1, &ram);
    vec_dst(s->g2, r->g2, &rbn);
    vec_dst(s->g1, r->g2, &ran);

    return 1 / rbm - 1 / ram - 1 / rbn + 1 / ran;
}

// Вычисление разности потенциалов на приемнике
double rec_dif(struct rec *r, struct vec *i)
{
    double d = 0;

    for (int s = 0; s < r->k.n; ++s) {
        d += r->k.dat[s] *
             i->dat[s] /
             (2 * M_PI * COND);
    }

    return d;
}

int main(int argc, char **argv)
{
    (void)argc;
    (void)argv;

    int sn = sizeof(sup) / sizeof(void *);
    int rn = sizeof(rec) / sizeof(void *);

    // Расчет коэффициентов

    for (int i = 0; i < rn; ++i) {
        if (vec_new(&rec[i]->k, sn)) {
            return -1;
        }

        for (int j = 0; j < sn; ++j) {
            rec[i]->k.dat[j] =
                k(&rec[i]->l, &sup[j]->l);
        }
    }

    struct vec in; // вектор сил токов
    struct vec id; // вектор приращения
    struct vec dc; // вектор измерений

    if (vec_new(&in, sn)) {
        return -1;
    }

    if (vec_new(&id, sn)) {
        return -1;
    }

    if (vec_new(&dc, rn)) {
        return -1;
    }

    in.dat[0] = 1.0;
    in.dat[1] = 1.2;
    in.dat[2] = 1.4;

    // Расчет синтетических данных (прямая задача)

    for (int i = 0; i < rn; ++i) {
        dc.dat[i] = rec_dif(rec[i], &in);
    }

    printf("d1 = %.7e\nd2 = %.7e\nd3 = %.7e\n",
        dc.dat[0], dc.dat[1], dc.dat[2]);

    struct imtx a; // матрица А
    struct imtx at;
    struct vec  b; // вектор правой части b
    struct vec  bt;

    if (mtx_new(&a, ((struct imtx_pps){
                        .n = sn, .m = sn}))) {
        return -1;
    }

    if (mtx_new(&at, ((struct imtx_pps){
                         .n = sn, .m = sn}))) {
        return -1;
    }

    if (vec_new(&b, sn)) {
        return -1;
    }

    if (vec_new(&bt, sn)) {
        return -1;
    }

    // Начальный вектор

    in.dat[0] = 0;
    in.dat[1] = 0;
    in.dat[2] = 0;

    for (int n = 0; n < 20; ++n) {
        // Вычисление функционала

        double f = 0;

        for (int k = 0; k < rn; ++k) {
            double dk = rec_dif(rec[k], &in);
            double ek = dk - dc.dat[k];
            double om = 1 / dc.dat[k];

            f += om * om * ek * ek;
        }

        printf("n = %d\nf = %.7e\ni1 = %.7e\n"
               "i2 = %.7e\ni3 = %.7e\n",
            n, f, in.dat[0], in.dat[1],
            in.dat[2]);

        if (f < 1e-17) {
            break;
        }

        for (int i = 0; i < sn; ++i) {
            for (int j = 0; j < sn; ++j) {
                double aij = 0;

                for (int k = 0; k < rn; ++k) {
                    double ki = rec[k]->k.dat[i];
                    double kj = rec[k]->k.dat[j];
                    double om = 1 / dc.dat[k];

                    aij += om *
                           om *
                           ki *
                           kj /
                           (4 *
                               M_PI *
                               M_PI *
                               COND *
                               COND);
                }

                a.dat[i][j] = aij;
            }

            double bi = 0;

            for (int k = 0; k < rn; ++k) {
                double ki = rec[k]->k.dat[i];
                double dk = rec_dif(rec[k], &in);
                double ek = dk - dc.dat[k];
                double om = 1 / dc.dat[k];

                bi += om *
                      om *
                      ek *
                      ki /
                      (2 * M_PI * COND);
            }

            b.dat[i] = -bi;
        }

        mtx_vdup(&a, &at);
        vec_dup(&b, &bt);

        if (dss_red_slv(&a, &id, &b)) {
            // Вырожденная матрица
            // Регуляризация

            double reg = 1e-10;

            for (int i = 0; i < sn; ++i) {
                at.dat[i][i] += reg;
            }

            if (dss_red_slv(&at, &id, &bt)) {
                goto end;
            }
        }

        // vec_cmb - комбинирование векторов

        vec_cmb(&in, &id, &in, 1);
    }

end:
    vec_cls(&in);
    vec_cls(&id);
    vec_cls(&dc);
    mtx_cls(&a);
    vec_cls(&b);

    for (int i = 0; i < rn; ++i) {
        vec_cls(&rec[i]->k);
    }

    return 0;
}
