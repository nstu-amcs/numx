#include <munit.h>

MunitResult test(const MunitParameter[], void *);

static char *c1_slv[] = {"fem", NULL};
static char *c1_msh[] = {"m1", NULL};
static char *c1_frq[] = {"f1", NULL};
static char *c1_usr[] = {"u1", NULL};
static char *c1_bdf[] = {"2", NULL};
static char *c1_ini[] = {"1", NULL};

static char *c2_slv[] = {"fem", NULL};
static char *c2_msh[] = {"m3", NULL};
static char *c2_frq[] = {"f1", NULL};
static char *c2_usr[] = {"u1", NULL};
static char *c2_bdf[] = {"2", NULL};
static char *c2_ini[] = {"1", NULL};

static MunitParameterEnum c1_pps[] = {
    {"slv", c1_slv},
    {"msh", c1_msh},
    {"frq", c1_frq},
    {"usr", c1_usr},
    {"bdf", c1_bdf},
    {"ini", c1_ini},
    {NULL,  NULL  },
};

static MunitParameterEnum c2_pps[] = {
    {"slv", c2_slv},
    {"msh", c2_msh},
    {"frq", c2_frq},
    {"usr", c2_usr},
    {"bdf", c2_bdf},
    {"ini", c2_ini},
    {NULL,  NULL  },
};

static void *c1_sup(const MunitParameter[], void *)
{
    return "pde/test/hmc/c1";
}

static void *c2_sup(const MunitParameter[], void *)
{
    return "pde/test/hmc/c2";
}

static MunitTest tests[] = {
    {"/c1", test, c1_sup, NULL, MUNIT_TEST_OPTION_NONE, c1_pps},
    {"/c2", test, c2_sup, NULL, MUNIT_TEST_OPTION_NONE, c2_pps},
    {NULL,  NULL, NULL,   NULL, MUNIT_TEST_OPTION_NONE, NULL  },
};

MunitSuite hmc_suite = {
    "/hmc",
    tests,                   // tests
    NULL,                    // suites
    1,                       // iterations
    MUNIT_SUITE_OPTION_NONE, // options
};

MunitSuite unit_pde_hmc_suite()
{
    return hmc_suite;
};
