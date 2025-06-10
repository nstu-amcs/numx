#include <munit.h>
#include <numx/vec/iss.h>

static MunitResult test_gmr(const MunitParameter[], void *)
{
    struct smtx m;
    struct vec  f;
    struct vec  x;

    mtx_new(&m, ((struct smtx_pps){.n = 3, .z = 2}));
    vec_new(&f, 3);
    vec_new(&x, 3);

    m.dr[0] = 1;
    m.dr[1] = 14;
    m.dr[2] = 20;
    m.lr[0] = -7;
    m.lr[1] = 0.1;
    m.ur[0] = 5;
    m.ur[1] = 8;

    m.ja[0] = 0;
    m.ja[1] = 0;

    m.ia[0] = 0;
    m.ia[1] = 0;
    m.ia[2] = 1;
    m.ia[3] = 2;

    f.dat[0] = 1;
    f.dat[1] = 2;
    f.dat[2] = 3;

    iss_gmr_slv(&m, &x, &f,
        &((struct iss_gmr_ops){
            .ops = {.err = 1e-10, .max = 20}
    }));

    mtx_cls(&m);
    vec_cls(&f);
    vec_cls(&x);

    return MUNIT_OK;
}

static MunitTest tests[] = {
    {"/gmr", test_gmr, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
    {NULL,   NULL,     NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
};

int main(int argc, char **argv)
{
    const MunitSuite suite = {
        "/numx/vec/iss/siss",
        tests,                   // tests
        NULL,                    // suites
        1,                       // iterations
        MUNIT_SUITE_OPTION_NONE, // options
    };

    return munit_suite_main(&suite, NULL, argc, argv);
}
