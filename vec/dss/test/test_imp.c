#include <errno.h>
#include <munit.h>
#include <stdio.h>
#include <string.h>
#include <numx/vec/dss.h>

static MunitResult test_red(const MunitParameter p[], void* c) {
  (void)p;
  (void)c;

  struct imtx m;
  struct vec x;
  struct vec f;

  struct imtx_pps pps = {3, 3};

  mtx_new(&m, pps);
  vec_new(&x, 3);
  vec_new(&f, 3);

  m.dat[0][0] = 4;
  m.dat[0][1] = 3;
  m.dat[0][2] = 5;
  m.dat[1][0] = 0;
  m.dat[1][1] = 0;
  m.dat[1][2] = 1;
  m.dat[2][0] = 2;
  m.dat[2][1] = -1;
  m.dat[2][2] = 2;

  f.dat[0] = 25;
  f.dat[1] = 3;
  f.dat[2] = 6;

  if (dss_red_slv(&m, &x, &f))
    fprintf(stderr, "err: %s\n", strerror(errno));

  munit_assert_double_equal(1.0, x.dat[0], 3);
  munit_assert_double_equal(2.0, x.dat[1], 3);
  munit_assert_double_equal(3.0, x.dat[2], 3);

  mtx_cls(&m);
  vec_cls(&x);
  vec_cls(&f);

  return MUNIT_OK;
}

static MunitTest itests[] = {
  {"/red", test_red, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  {NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
};

static const MunitSuite suite = {"/imp", itests, NULL, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char** argv) {
  return munit_suite_main(&suite, NULL, argc, argv);
}
