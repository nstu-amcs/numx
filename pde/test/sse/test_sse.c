#include <munit.h>
#include <numx/pde/sse.h>
#include <stdio.h>

// clang-format off

static void* test_fde_sup(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  return NULL;
}

static void test_fde_tdn(void* c) {
  (void)c;
}

static MunitResult test_fde(const MunitParameter pps[], void* ctx) { 
  (void)pps;
  (void)ctx;

  return MUNIT_OK; 
}

static char* cases[] = { "case1", NULL };
static MunitParameterEnum pms[] = {{ "case", cases }, { NULL, NULL }};

static MunitTest tests[] = {
  { "/fde", test_fde, test_fde_sup, test_fde_tdn, MUNIT_TEST_OPTION_NONE, pms},
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};

static const MunitSuite suite = {
  "/pde/sse", tests, NULL, 1, MUNIT_SUITE_OPTION_NONE,
};

int main(int argc, char** argv) {
  return munit_suite_main(&suite, NULL, argc, argv);
}
