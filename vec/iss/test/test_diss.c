#include <munit.h>
#include <numx/vec/iss.h>
#include <stdarg.h>
#include <stdio.h>

// clang-format off

//     | 1 2 0 0 1 |
//     | 0 7 3 0 0 |
// m = | 0 0 9 3 0 |, x = (1, 2, 3, 4, 5)
//     | 1 0 0 4 2 |
//     | 0 5 0 0 8 |
static MunitResult test_diss_rlx(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  struct dmtx m;
  struct vec x;
  struct vec f;

  munit_assert_int(0, ==, mtx_new(&m, ((struct dmtx_pps){.n = 5, .d = 4})));

  m.la[0] = 0;
  m.la[1] = 1;
  m.la[2] = 4;
  m.la[3] = -3;

  m.ad[0][0] = 1; m.ad[0][1] = 2; m.ad[0][2] = 1; m.ad[0][3] = 0;
  m.ad[1][0] = 7; m.ad[1][1] = 3; m.ad[1][2] = 0; m.ad[1][3] = 0;
  m.ad[2][0] = 9; m.ad[2][1] = 3; m.ad[2][2] = 0; m.ad[2][3] = 0;
  m.ad[3][0] = 4; m.ad[3][1] = 2; m.ad[3][2] = 0; m.ad[3][3] = 1;
  m.ad[4][0] = 8; m.ad[4][1] = 0; m.ad[4][2] = 0; m.ad[4][3] = 5;

  munit_assert_int(0, ==, vec_new(&x, 5, 1.0, 2.0, 3.0, 4.0, 5.0));
  munit_assert_int(0, ==, vec_new(&f, 5));

  struct iss_rlx_ops ops = {
      .ops =
          {
              .eps = 1e-15,
              .max = 1000,
          },
      .omg = 1.2,
  };

  struct iss_itr i;

  ops.ops.itr.ctx = &i;
  ops.ops.itr.call = &iss_itr_cap;

  munit_assert_int(0, ==, mtx_vmlt(&m, &x, &f));
  munit_assert_int(0, ==, vec_rst(&x));
  munit_assert_int(0, ==, iss_rlx_slv(&m, &x, &f, ops));

  munit_assert_double_equal(1.0, x.dat[0], 12);
  munit_assert_double_equal(2.0, x.dat[1], 12);
  munit_assert_double_equal(3.0, x.dat[2], 12);
  munit_assert_double_equal(4.0, x.dat[3], 12);
  munit_assert_double_equal(5.0, x.dat[4], 12);

  printf("[rlx] num: %d; err: %.7e\n", i.num, i.err);

  mtx_cls(&m);
  vec_cls(&x);
  vec_cls(&f);

  return MUNIT_OK;
}

//     | 1 2 0 0 1 |
//     | 0 7 3 0 0 |
// m = | 0 0 9 3 0 |, x = (1, 2, 3, 4, 5)
//     | 1 0 0 4 2 |
//     | 0 5 0 0 8 |
static MunitResult test_diss_jac(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  struct dmtx m;
  struct vec x;
  struct vec f;

  munit_assert_int(0, ==, mtx_new(&m, ((struct dmtx_pps){.n = 5, .d = 4})));

  m.la[0] = 0;
  m.la[1] = 1;
  m.la[2] = 4;
  m.la[3] = -3;

  m.ad[0][0] = 1; m.ad[0][1] = 2; m.ad[0][2] = 1; m.ad[0][3] = 0;
  m.ad[1][0] = 7; m.ad[1][1] = 3; m.ad[1][2] = 0; m.ad[1][3] = 0;
  m.ad[2][0] = 9; m.ad[2][1] = 3; m.ad[2][2] = 0; m.ad[2][3] = 0;
  m.ad[3][0] = 4; m.ad[3][1] = 2; m.ad[3][2] = 0; m.ad[3][3] = 1;
  m.ad[4][0] = 8; m.ad[4][1] = 0; m.ad[4][2] = 0; m.ad[4][3] = 5;

  munit_assert_int(0, ==, vec_new(&x, 5, 1.0, 2.0, 3.0, 4.0, 5.0));
  munit_assert_int(0, ==, vec_new(&f, 5));

  struct iss_jac_ops ops = {
      .ops =
          {
              .eps = 1e-15,
              .max = 1000,
          },
      .omg = 1.2,
  };

  struct iss_itr i;

  ops.ops.itr.ctx = &i;
  ops.ops.itr.call = &iss_itr_cap;

  munit_assert_int(0, ==, mtx_vmlt(&m, &x, &f));
  munit_assert_int(0, ==, vec_rst(&x));
  munit_assert_int(0, ==, iss_jac_slv(&m, &x, &f, ops));

  munit_assert_double_equal(1.0, x.dat[0], 12);
  munit_assert_double_equal(2.0, x.dat[1], 12);
  munit_assert_double_equal(3.0, x.dat[2], 12);
  munit_assert_double_equal(4.0, x.dat[3], 12);
  munit_assert_double_equal(5.0, x.dat[4], 12);

  printf("[rlx] num: %d; err: %.7e\n", i.num, i.err);

  mtx_cls(&m);
  vec_cls(&x);
  vec_cls(&f);

  return MUNIT_OK;
}

static MunitTest tests[] = {
  {"/rlx", test_diss_rlx, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  {"/jac", test_diss_jac, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  // {"/sor", test_diss_sor_all, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  {NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL}
};

static MunitSuite suites[] = {
    {"/diss", tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
    {NULL, NULL, NULL, 1, MUNIT_SUITE_OPTION_NONE},
};

static const MunitSuite suite = {"vec/iss", NULL, suites, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char** argv) { return munit_suite_main(&suite, NULL, argc, argv); }