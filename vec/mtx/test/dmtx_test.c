#include <munit.h>
#include <numx/vec/mtx.h>

//     | 1 2 0 |
// m = | 0 7 3 |, x = (1, 2, 3), mx = f = (5, 23, 0)
//     | 0 0 0 |
static MunitResult test_dmtx_vmlt_asym(const MunitParameter[], void*) {
  struct dmtx m;
  struct vec x;
  struct vec f;

  munit_assert_int(0, ==, mtx_new(&m, ((struct dmtx_pps){3, 2})));

  m.la[0] = 0;
  m.la[1] = 1;

  m.ad[0][0] = 1;
  m.ad[1][0] = 7;
  m.ad[2][0] = 0;
  m.ad[0][1] = 2;
  m.ad[1][1] = 3;
  m.ad[2][1] = 0;

  munit_assert_int(0, ==, vec_new(&x, 3, 1.0, 2.0, 3.0));
  munit_assert_int(0, ==, vec_new(&f, 3));

  munit_assert_int(0, ==, mtx_vmlt(&m, &x, &f));

  munit_assert_double_equal(5.0, f.dat[0], 15);
  munit_assert_double_equal(23.0, f.dat[1], 15);
  munit_assert_double_equal(0.0, f.dat[2], 15);

  mtx_cls(&m);
  vec_cls(&x);
  vec_cls(&f);

  return MUNIT_OK;
}

static MunitTest vmlt_tests[] = {{"/asym", test_dmtx_vmlt_asym, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
                                 {NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL}};

static MunitSuite suites[] = {
    {"/vmlt", vmlt_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
    {NULL, NULL, NULL, 1, MUNIT_SUITE_OPTION_NONE},
};

static const MunitSuite suite = {"vec/mtx/dmtx", NULL, suites, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char** argv) { return munit_suite_main(&suite, NULL, argc, argv); }
