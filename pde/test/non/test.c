#include <munit.h>

static MunitResult test_fem(const MunitParameter[], void *)
{
    return MUNIT_OK;
}

static char *msh[] = {"m1", "m2", NULL};
static char *usr[] = {"u1", "u2", NULL};

static MunitParameterEnum pps[] = {
    {"msh", msh },
    {"usr", usr },
    {NULL,  NULL},
};

static MunitTest tests[] = {
    {"/fem", test_fem, NULL, NULL, MUNIT_TEST_OPTION_NONE, pps}
};

MunitSuite pbc_suite = {
    "/pbc",
    tests,                   // tests
    NULL,                    // suites
    1,                       // iterations
    MUNIT_SUITE_OPTION_NONE, // options
};

MunitSuite unit_pde_pbc_suite()
{
    return pbc_suite;
}
