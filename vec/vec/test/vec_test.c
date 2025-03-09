#include <math.h>
#include <munit.h>
#include <numx/vec/vec.h>
#include <stdio.h>
#include <string.h>

static MunitResult test_new(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  struct vec v;

  vec_new_ini(&v, 3, 1.0, 2.0, 3.0);

  munit_assert_double_equal(1.0, v.dat[0], 15);
  munit_assert_double_equal(2.0, v.dat[1], 15);
  munit_assert_double_equal(3.0, v.dat[2], 15);

  vec_cls(&v);

  return MUNIT_OK;
}

static MunitResult test_cmb(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  struct vec v1;
  struct vec v2;

  vec_new_ini(&v1, 3, 1.0, 2.0, 3.0);
  vec_new_ini(&v2, 3, 9.0, -1.0, 20.0);

  munit_assert_int(0, ==, vec_cmb(&v2, &v1, &v1, -1));

  munit_assert_double_equal(8.0, v1.dat[0], 15);
  munit_assert_double_equal(-3.0, v1.dat[1], 15);
  munit_assert_double_equal(17.0, v1.dat[2], 15);

  vec_cls(&v1);
  vec_cls(&v2);

  return MUNIT_OK;
}

static MunitResult test_nrm(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  struct vec v;
  double nrm = 0;

  vec_new_ini(&v, 3, 1.0, 2.0, 3.0);

  munit_assert_int(0, ==, vec_nrm(&v, &nrm));
  munit_assert_double_equal(sqrt(14.0), nrm, 15);

  vec_cls(&v);

  return MUNIT_OK;
}

static MunitTest tests[] = {
    {"/new", test_new, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
    {"/cmb", test_cmb, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
    {"/nrm", test_nrm, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
    {NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
};

static const MunitSuite suite = {"/vec/vec", tests, NULL, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char** argv) { return munit_suite_main(&suite, NULL, argc, argv); }
