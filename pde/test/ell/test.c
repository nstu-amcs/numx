#include <munit.h>

static MunitResult test_fem(const MunitParameter[], void *)
{
    return MUNIT_OK;
}

static MunitTest tests[] = {
    {"/fem", test_fem, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL}
};

MunitSuite ell_suite = {
    "/ell",
    tests,                   // tests
    NULL,                    // suites
    1,                       // iterations
    MUNIT_SUITE_OPTION_NONE, // options
};

MunitSuite unit_pde_ell_suite()
{
    return ell_suite;
}
