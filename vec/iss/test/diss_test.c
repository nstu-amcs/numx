#include <munit.h>
#include <numx/vec/iss.h>
#include <stdarg.h>

void itr_log(void*, int n, ...) {
  va_list arg;
  va_start(arg, n);

  if (n != 2)
    return;

  int i = va_arg(arg, int);
  double r = va_arg(arg, double);

  munit_logf(MUNIT_LOG_INFO, "itr: %d, res: %.5e", i, r);

  va_end(arg);
}

//     | 1 2 0 |
// m = | 0 7 3 |, x = (1, 2, 3), mx = f = (5, 23, 27)
//     | 0 0 9 |
static MunitResult test_diss_rlx_posd(const MunitParameter[], void*) {
  struct dmtx m;
  struct vec x;
  struct vec f;

  munit_assert_int(0, ==, mtx_new(&m, ((struct dmtx_pps){3, 2})));

  m.la[0] = 0;
  m.la[1] = 1;

  m.ad[0][0] = 1;
  m.ad[1][0] = 7;
  m.ad[2][0] = 9;
  m.ad[0][1] = 2;
  m.ad[1][1] = 3;
  m.ad[2][1] = 0;

  munit_assert_int(0, ==, vec_new(&x, 3, 1.0, 2.0, 3.0));
  munit_assert_int(0, ==, vec_new(&f, 3));

  struct iss_rlx_opt opt = {
      .opt =
          {
              .eps = 1e-10,
              .max = 1000,
          },
      .omg = 1,
  };

  opt.opt.itr.call = &itr_log;

  munit_assert_int(0, ==, mtx_vmlt(&m, &x, &f));
  munit_assert_int(0, ==, vec_rst(&x));
  munit_assert_int(0, ==, iss_rlx_slv(&m, &x, &f, opt));

  munit_assert_double_equal(1.0, x.dat[0], 5);
  munit_assert_double_equal(2.0, x.dat[1], 5);
  munit_assert_double_equal(3.0, x.dat[2], 5);

  mtx_cls(&m);
  vec_cls(&x);
  vec_cls(&f);

  return MUNIT_OK;
}

static MunitTest vmlt_tests[] = {{"/posd", test_diss_rlx_posd, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
                                 {NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL}};

static MunitSuite suites[] = {
    {"/rlx", vmlt_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
    {NULL, NULL, NULL, 1, MUNIT_SUITE_OPTION_NONE},
};

static const MunitSuite suite = {"vec/iss/diss", NULL, suites, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char** argv) { return munit_suite_main(&suite, NULL, argc, argv); }