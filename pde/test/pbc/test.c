#include <munit.h>

MunitResult test(const MunitParameter[], void *);

static char *c1_slv[] = {"fem", NULL};
static char *c1_msh[] = {"m1", NULL};
static char *c1_frq[] = {"f1", NULL};
static char *c1_usr[] = {"u1", NULL};
static char *c1_bdf[] = {"4", NULL};
static char *c1_ini[] = {"3", NULL};

static MunitParameterEnum c1_pps[] = {
    {"slv", c1_slv},
    {"msh", c1_msh},
    {"frq", c1_frq},
    {"usr", c1_usr},
    {"bdf", c1_bdf},
    {"ini", c1_ini},
    {NULL,  NULL  },
};

static void *c1_sup(const MunitParameter[], void *)
{
    return "pde/test/pbc/c1";
}

static MunitTest tests[] = {
    {"/c1", test, c1_sup, NULL, MUNIT_TEST_OPTION_NONE, c1_pps},
    {NULL,  NULL, NULL,   NULL, MUNIT_TEST_OPTION_NONE, NULL  },
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
};
